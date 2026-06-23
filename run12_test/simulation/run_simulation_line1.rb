#! /usr/bin/env ruby

require 'comptonsoft'

def run_simulation(num, random, output)
  energy = 1332.5

  sim = ComptonSoft::Simulation.new
  sim.output = output
  sim.random_seed = random
  sim.verbose = 0
  sim.print_detector_info
  sim.set_database(detector_configuration: "../database/detector_configuration.xml",
                   detector_parameters: "../database/detector_parameters.xml")
  sim.set_gdml "../database/mass_model.gdml"

  sim.set_physics(hadron_hp: false, cut_value: 0.001)
  sim.enable_timing_process

  sim.set_physics(hadron_hp: false, cut_value: 0.001, radioactive_decay: true)
  #sim.set_primary_generator :NucleusPrimaryGen, {
  #  atomic_number: 27,
  #  mass_number: 60,
  #  energy: 0.0,
  #  position: vec(159.4, 0.0, 37.5),
  #}
  sim.set_primary_generator :PointSourcePrimaryGen, {
    particle: "gamma",
    spectral_distribution: "gaussian",
    energy_mean: energy,
    energy_sigma: energy*0.001,
    # position: vec(14.5, 0.0, -12.25),
    position: vec(159.4, 0.0, 37.5),
    direction: vec(-1.0, 0.0, 0.0),
    theta_min: 0.0,
    theta_max: 90.0*Math::PI/180.0,
  }

  sim.visualize(mode: 'OGLSQt')
  sim.run(num)
end

### main ###
#num = 1000000000
num = 10
random = 0
output = "simulation_line1.root"
run_simulation(num, random, output)
