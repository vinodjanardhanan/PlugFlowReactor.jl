module PlugFlowReactor



using LightXML, Printf, DifferentialEquations, Sundials
using SurfaceReactions, GasphaseReactions, ReactionCommons, IdealGas, RxnHelperUtils

export plug

global os_streams


"""
This is the calling function for executing the plug reactor with userdefined rate calculation
#   Method-1
plug(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
-   input_file: the xml input file for batch reactor
-   lib_dir: the direcrtory in which the data files are present. It must be the relative path
-   user_defined: A function which calculates the species source terms, to be supplied by the user  
"""
function plug(input_file::AbstractString, lib_dir::AbstractString, user_defined::Function; sens= false)    
    chem = Chemistry(false, false, true, user_defined)
    return plug(input_file, lib_dir, sens, chem)
end


"""
This is the calling function for executing the plug reactor with chemistry input files 
#  Method-2
plug(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
-   input_file: the xml input file plug flow reactor configuration
-   lib_dir: the direcrtory in which the data files are present (mechanism file and therm.dat file). 
-   sens: sensitivity calculation 
-   surfchem : surface chemistry; false implies no surface chemistry calculation 
-   gaschem : gasphase chemistry; false implies no gasphase chemistry calculation 
"""
function plug(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
    chem = Chemistry(surfchem, gaschem, false, f->())
    return plug(input_file, lib_dir, sens, chem)
    
end



#=
This is the general iterface for starting the plug flow code by reading the input files.   
=#
function plug(input_file::AbstractString, lib_dir::AbstractString, sens, chem::Chemistry)
    
    #locate the lib directory    
    thermo_file = "therm.dat"
    thermo_file = get_path(lib_dir, thermo_file)


    local gasphase = Array{String,1}

    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)

    """ In case the gasphase mechanism is present then the gasphase species are read from the 
     gasphase mechanism file and not from the xml input. So this must be the first action before any
     other input parameters are read 
     """
    if chem.gaschem
        mech_file = get_text_from_xml(xmlroot,"gas_mech")
        mech_file = get_path(lib_dir, mech_file)
        gmd = compile_gaschemistry(mech_file)        
        gasphase = gmd.gm.species
    end


    #Get gasphase species list from the xml input
    if !chem.gaschem
        gasphase = get_collection_from_xml(xmlroot,"gasphase")
    end
    thermo_obj = IdealGas.create_thermo(gasphase,thermo_file)        

    #Get the molefractions
    mole_fracs = get_molefraction_from_xml(xmlroot, thermo_obj.molwt, gasphase)
    
    
    #Convert mole fractions to mass fractions
    mass_fracs = zeros( length(gasphase))
    molefrac_to_massfrac!(mass_fracs,mole_fracs,thermo_obj.molwt)
    

    #Read the inlet temperature
    local T = get_value_from_xml(xmlroot,"T")

    #Read the velocity 
    local u = get_value_from_xml(xmlroot,"u")

    #Read the pressure
    local p = get_value_from_xml(xmlroot,"p")

    #Read the length of the reactor
    local l = get_value_from_xml(xmlroot,"length")

    #Read the diameter of the reactor 
    local dia = get_value_from_xml(xmlroot,"dia")
    local Ac = π*dia^2/4.0
    local As_per_unit_length = π*dia

    #Read the wall temperature
    local Tw = get_value_from_xml(xmlroot,"Tw")

    #Read the surface area per unit length
    local cat_geom = get_value_from_xml(xmlroot,"cat-geom-factor")

    #Read the temperature condition     
    isothermal = lowercase(get_text_from_xml(xmlroot,"isothermal")) == "true" ? true : false
    
    #Create the mechanism definition    
    
    if chem.surfchem
        #Get the mechanism file from xml
        mech_file = get_text_from_xml(xmlroot,"surface_mech")
        # mech_file = lib_dir*"/"*mech_file
        mech_file = get_path(lib_dir, mech_file)
        md = SurfaceReactions.compile_mech(mech_file,thermo_obj,gasphase)
    end
        
    

    #create the soln vector
    soln = deepcopy(mass_fracs) #mass fractions
    push!(soln,density(mole_fracs,thermo_obj.molwt,T,p)*u) #mass flux ρu
    if !isothermal
        push!(soln,T)
    end

    #file output streams for saving the data
    g_stream = open(output_file(input_file, "gas_profile.dat"),"w")
    s_stream = open(output_file(input_file, "surf_covg.dat"),"w")
    csv_g_stream = open(output_file(input_file, "gas_profile.csv"),"w")
    csv_s_stream = open(output_file(input_file, "surf_covg.csv"),"w")
    global os_streams = (g_stream,s_stream, csv_g_stream, csv_s_stream)
    geometry = (Ac,As_per_unit_length,cat_geom)    

    create_header(g_stream,["z","T","p","u","rho"],gasphase)
    write_csv(csv_g_stream,["z","T","p","u","rho"],gasphase)
    if chem.surfchem
        create_header(s_stream,"z","T",md.sm.species)
        write_csv(csv_s_stream,"z","T",md.sm.species)
    end
    #Set up the problem
    t_span = (0,l)
    n_species = length(gasphase)
    if chem.surfchem
        covg = md.sm.si.ini_covg
        n_species  += length(md.sm.species)
    end
    
    rate = zeros(n_species)      
    all_conc = zeros(n_species)      

    # create the reaction state object for the calculation of surface reaction rates 
    if chem.surfchem
        #Create the state object    
        surf_conc = similar(covg)
        rxn_rate = zeros(length(md.sm.reactions))
        sr_state = SurfaceRxnState(T, p, mole_fracs, covg, surf_conc ,  rxn_rate , rate, all_conc)
        #this step is to get the steady state coverage before starting the plug flow integration
        t, sr_state = SurfaceReactions.calculate_ss_molar_production_rates!(sr_state,thermo_obj,md,1.0)
        # params = (state, thermo_obj, md, geometry, os_streams,chem)
    end

    # create the reaction state for gasphase 
    if chem.gaschem
        conc = similar(mole_fracs)
        source = zeros(length(conc))
        g_all = similar(mole_fracs)    
        Kp = zeros(length(gmd.gm.reactions))            
        rxn_rate = zeros(length(Kp))        
        gs_state = GasphaseState(T, p, mole_fracs, conc, rxn_rate, source, g_all, Kp)
        
    end

    if chem.userchem && !chem.surfchem && !chem.gaschem        
        source = zeros(length(mole_fracs))
        usr_state = UserDefinedState(T,p, mole_fracs, thermo_obj.molwt,gasphase, source)
    end


    if chem.surfchem && !chem.gaschem
        params = (s_state = sr_state, thermo = thermo_obj, smd = md, geometry= geometry, chem = chem)
    elseif chem.gaschem && !chem.surfchem
        params = (g_state = gs_state, thermo = thermo_obj, gmd = gmd, geometry = geometry,  chem=chem)        
    elseif chem.gaschem && chem.surfchem
        params = (g_state = gs_state, s_state=sr_state, thermo =thermo_obj, smd=md, gmd=gmd, geometry=geometry, chem=chem)
    elseif chem.userchem
        params = (state = usr_state, thermo = thermo_obj, md = nothing, geometry = geometry, chem = chem)
    end


    prob = ODEProblem(residual!,soln,t_span,params)
    if sens == true
        return (params, prob, t_span)
    end

    alg = CVODE_BDF()    
    cb = FunctionCallingCallback(write_out_put)    
    sol = solve(prob, alg, reltol=1e-6, abstol=1e-10, save_everystep=false,callback=cb);        
    
    
    close(s_stream)
    close(g_stream)    
    close(csv_g_stream)
    close(csv_s_stream)
    return Symbol(sol.retcode)
end


"""
A function for call from other packages, mainly intended for reactor network modeling 
#   Method-3
plug(inlet_comp, T, p, vel, l; chem, thermo_obj, md, geometry)
-   inlet_comp : Dictionary of species name and mole fractions. In the case of gasphase chemistry the inlet_comp 
may only contain species that are present at the inlet 
-   T : Temperature in K
-   p : Pressure in Pa 
-   vel : inlet velocity in m/s 
-   l : length of the reactor in m
-   chem: Chemistry to be invoked: tuple Chemistry(surfchem, , true, user_defined)
-   thermo_obj : Thermo Object 
-   md : Mechanism definition(Please refer SurfaceReactions.jl and Ga
-   geometry : tuple of reactor geometry (dia,cat_geom)      
    -   dia: reactor diameter 
    -   cat_geom : enhancement factor for the reaction rate due to the catalyst geometry effect
"""
function plug(inlet_comp, T, p, vel, l; chem, thermo_obj, md, geometry)

    dia = geometry[Symbol("dia")]    
    Ac = π*dia^2/4.0
    As_per_unit_length = π*dia
    cat_geom = geometry[Symbol("cat_geom")]
    geom = (Ac,As_per_unit_length,cat_geom)    

    function  get_mole_fracs(species, inlet_comp)
        mole_fracs = zeros(length(species))
        for k in eachindex(species)
            if in(species[k], collect(keys(inlet_comp)))
                mole_fracs[k] = inlet_comp[species[k]]
            end
        end 
        mole_fracs       
    end

    if chem.surfchem
        species = collect(keys(inlet_comp))
        mole_fracs = get_mole_fracs(species, inlet_comp)
        covg = md.sm.si.ini_covg
        #Create the state object    
        surf_conc = similar(covg)
        rxn_rate = zeros(length(md.sm.reactions))
        n_species = length(mole_fracs) + length(covg)
        rate = zeros(n_species)
        all_conc = zeros(n_species)
        sr_state = SurfaceRxnState(T, p, mole_fracs, covg, surf_conc, rxn_rate, rate, all_conc)
        #this step is to get the steady state coverage before starting the plug flow integration
        t, sr_state = SurfaceReactions.calculate_ss_molar_production_rates!(sr_state,thermo_obj,md,1.0)   
            
        params = (s_state = sr_state, thermo = thermo_obj, smd = md, geometry = geom, chem = chem)
    end

    if chem.gaschem        
        species = md.gm.species
        mole_fracs = get_mole_fracs(md.gm.species, inlet_comp)
        conc = zeros(length(mole_fracs))
        source = zeros(length(conc))
        g_all = zeros(length(mole_fracs))
        Kp = zeros(length(md.gm.reactions))            
        rxn_rate = zeros(length(Kp))                

        gs_state = GasphaseState(T, p, mole_fracs, conc, rxn_rate, source, g_all, Kp)       
        params = (g_state = gs_state, thermo = thermo_obj, gmd = md, geometry= geom, chem = chem)              
    end


    ρ = density(mole_fracs,thermo_obj.molwt,T,p)
    soln = similar(mole_fracs)
    molefrac_to_massfrac!(soln, mole_fracs, thermo_obj.molwt)    
    push!(soln, vel*ρ)

    
    t_span = (0, l)
    prob = ODEProblem(residual!, soln, t_span, params)
    alg = CVODE_BDF()
    sol = solve(prob, alg , reltol = 1e-6, abstol=1e-10, save_everystep=false)
    mass_fracs = sol.u[end][1:end-1]
    mole_fracs = zeros(length(mass_fracs))
    massfrac_to_molefrac!(mole_fracs, mass_fracs, thermo_obj.molwt)    
    return sol.t, Dict(species .=> mole_fracs)    
end



#=
residual function that defined the governing equations    
=#
function residual!(du,u,p,t)    

    # state,thermo_obj,md,geom,os_stream = p
    geom = p[:geometry]
    Ac,Aspul,cat_geom = geom
    
    #Convert massfractions to molefractions
    if p[:chem].gaschem && p[:chem].surfchem
        ng = length(p[:g_state].mole_frac)    
    elseif p[:chem].surfchem 
        ng = length(p[:s_state].mole_frac)    
    elseif p[:chem].gaschem
        ng = length(p[:g_state].mole_frac)    
    elseif p[:chem].userchem
        ng = length(p[:state].mole_frac)
    end

    # massfracs = u[1:ng]
    
    
    
    if p[:chem].surfchem 
        #Calculate the molar production rates due to surface reactions 
        massfrac_to_molefrac!(p[:s_state].mole_frac ,u[1:ng], p[:thermo].molwt)
        SurfaceReactions.calculate_ss_molar_production_rates!(p[:s_state],p[:thermo],p[:smd],1.0)
        p[:s_state].source[1:ng] *= cat_geom*(Aspul/Ac)
    end

    if p[:chem].gaschem && ! p[:chem].surfchem        
        massfrac_to_molefrac!(p[:g_state].mole_frac ,u[1:ng], p[:thermo].molwt)
        #Calculate the molar production rates due to gasphase reactions 
        GasphaseReactions.calculate_molar_production_rates!(p[:g_state], p[:gmd], p[:thermo])
    end

    if p[:chem].userchem
        p[:chem].udf(p[:state])
    end
    
    #species residual
    if p[:chem].surfchem && !p[:chem].gaschem            
        du[1:ng] = (p[:s_state].source[1:ng] .* p[:thermo].molwt)/u[ng+1]
        #mass flux
        du[ng+1] = sum(p[:s_state].source[1:ng] .* p[:thermo].molwt)
    elseif !p[:chem].surfchem && p[:chem].gaschem            
        du[1:ng] = (p[:g_state].source[1:ng] .* p[:thermo].molwt)/u[ng+1]
        #mass flux
        du[ng+1] = sum(p[:g_state].source[1:ng] .* p[:thermo].molwt)
    elseif p[:chem].surfchem && p[:chem].gaschem
        du[1:ng] = ((p[:g_state].source[1:ng] + p[:s_state].source[1:ng]) .* p[:thermo].molwt)/u[ng+1]
        #mass flux
        du[ng+1] = sum(p[:g_state].source[1:ng] .* p[:thermo].molwt) + sum(p[:s_state].source[1:ng] .* p[:thermo].molwt)
    elseif p[:chem].userchem
        du[1:ng] = p[:state].source[1:ng] .* p[:thermo].molwt/ u[ng+1]
    end

    
end




function write_out_put(u,t,integrator)    
    # state = integrator.p[1]
    # thermo_obj = integrator.p[2]
    # g_stream, s_stream = integrator.p[5]
    if integrator.p[:chem].surfchem && !integrator.p[:chem].gaschem
        state = integrator.p[:s_state]
    elseif integrator.p[:chem].gaschem
        state = integrator.p[:g_state]
    elseif integrator.p[:chem].userchem
        state = integrator.p[:state]
    end
    g_stream, s_stream, csv_g_stream, csv_s_stream = os_streams    
    d = density(state.mole_frac,integrator.p[:thermo].molwt,state.T,state.p)
    vel = u[length(state.mole_frac)+1]/d
    write_to_file(g_stream,t,state.T,state.p,vel,d,state.mole_frac)
    write_csv(csv_g_stream,t,state.T,state.p,vel,d,state.mole_frac)
    if integrator.p[:chem].surfchem
        write_to_file(s_stream,t,integrator.p[:s_state].T,integrator.p[:s_state].covg)
        write_csv(csv_s_stream,t,integrator.p[:s_state].T,integrator.p[:s_state].covg)
    end
    @printf("%.4e\n", t)   
end




end
