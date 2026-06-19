#! /usr/bin/env ruby

require 'comptonsoft'

def run_simulation(num, za_pair, random, output)

  sim = ComptonSoft::Simulation.new
  sim.random_seed = random
  sim.output      = output
  sim.verbose     = 0
  sim.print_detector_info

  sim.set_gdml "../database/mass_model_PbCollimator_run11_actual.gdml"
  sim.set_database(detector_configuration: "../database/detector_configuration.xml",
                   detector_parameters: "../database/detector_parameters.xml")
  sim.set_physics(radioactive_decay: true)
  sim.set_primary_generator :NucleusPrimaryGen, {
    atomic_number: za_pair[0],
    mass_number: za_pair[1],
    energy: 0.0,
    position: vec(32.0+10.0+14.5, 0.0, -12.5),
  }

  sim.run(num)
end

### main ###
dirOutput = "products"
num       = ARGV[0].to_i
random    = ARGV[1].to_i
za_pair   = [ARGV[2].to_i, ARGV[3].to_i]

output    = "#{dirOutput}/simulation_#{num}_Z#{za_pair[0]}_A#{za_pair[1]}_seed#{random}.root"
run_simulation(num, za_pair, random, output)
