#! /usr/bin/env ruby

require 'comptonsoft'

def run_simulation(num, energy, random)

  sim = ComptonSoft::Simulation.new
  sim.random_seed = random
  sim.verbose = 0
  sim.print_detector_info
  sim.set_database(detector_configuration: "../database/detector_configuration.xml",
                   detector_parameters: "../database/detector_parameters.xml")
  sim.set_gdml "../database/mass_model_PbCollimator_run11_actual.gdml"

  sim.set_physics(hadron_hp: false, cut_value: 0.001)
  sim.enable_timing_process

  sim.set_primary_generator :PointSourcePrimaryGen, {
    particle: "gamma",
    spectral_distribution: "gaussian",
    energy_mean: energy,
    energy_sigma: 0.0,
    position: vec(32.0+10.0+14.5, 0.0, -12.5),
    direction: vec(-1.0, 0.0, 0.0),
    theta_min: 0.0,
    theta_max: 90.0*Math::PI/180.0,
  }

  sim.visualize(mode: 'OGLSQt')
  sim.run(num)
end

### main ###
num = 10
#energy = 1332.5
energy = 100.0
random = 0
run_simulation(num, energy, random)
