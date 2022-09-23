using PlugFlowReactor
using Test
using IdealGas, RxnHelperUtils, SurfaceReactions, ReactionCommons, GasphaseReactions

@testset "PlugFlowReactor.jl" begin

    if Sys.isapple()  || Sys.islinux()
        lib_dir = "lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end

    @testset "Testing surface chemistry" begin        
        input_file = joinpath("plug_surf", "plug.xml")
        retcode = plug(input_file, lib_dir, surfchem=true)
        @test retcode == Symbol("Success")
    end

    @testset "Testing gas chemistry of h2 + o2" begin
        
        input_file = joinpath("plug_h2o2", "plug.xml")

        retcode = plug(input_file, lib_dir, gaschem=true)
        @test retcode == Symbol("Success")
    end

    @testset "Testing gas chemistry using grimech" begin
        input_file = joinpath("plug_ch4", "plug.xml")
        retcode = plug(input_file, lib_dir, gaschem=true)
        @test retcode == Symbol("Success")
    end

    @testset "Testing surface chemistry with interface call " begin
        inlet_comp = Dict("CH4"=>0.25,"H2O"=>0.0, "H2"=>0.0, "CO"=>0.0, "CO2"=>0.25, "O2"=>0.0, "N2"=>0.5)
        T = 1073.15
        p = 1e5
        u = 0.1
        l = 0.3
        gasphase = collect(keys(inlet_comp))
        thermo_obj = IdealGas.create_thermo(gasphase, get_path(lib_dir,"therm.dat"))
        md = SurfaceReactions.compile_mech(get_path(lib_dir,"ch4ni.xml"),thermo_obj,gasphase)
        chem = Chemistry(true, false, false, f->())
        geometry = (dia=0.005,cat_geom=1000.0)        
        retcodes = plug(inlet_comp, T, p, u, l; chem=chem, thermo_obj=thermo_obj, md=md, geometry=geometry)                       
        @test retcodes[1][end] == l
    end

    @testset "Testing gasphase chemistry with interface call " begin
        
        mech_file = get_path(lib_dir, "h2o2.dat")
        gmd = compile_gaschemistry(mech_file)        
        gasphase = gmd.gm.species
        thermo_obj = IdealGas.create_thermo(gasphase, get_path(lib_dir,"therm.dat"))        
        inlet_comp = Dict("O2" => 0.25,"N2" => 0.5, "H2" => 0.25)
        
        T = 1073.15
        p = 1e5
        u = 0.1
        l = 0.3
                
        chem = Chemistry(false, true, false, f->())      
        geometry = (dia=0.005,cat_geom=1.0)               
        retcodes = plug(inlet_comp, T, p, u, l; chem=chem, thermo_obj=thermo_obj, md=gmd, geometry=geometry)                                
        @test retcodes[1][end] == l
    end


    @testset "Testing user defined chemistry " begin        
        input_file = joinpath("plug_udf", "plug.xml")
        function udf(state)
            state.source[1:end] .= 0.0            
        end
        retcode = plug(input_file, lib_dir, udf)
        @test retcode == Symbol("Success")        
    end
end
