#! /usr/bin/env ruby

require 'comptonsoft'


def run_simulation(num, energy, detector_parameter_file, random, output)
  script_dir = File.expand_path(__dir__)

  sim = ComptonSoft::Simulation.new
  sim.output      = output
  sim.random_seed = random
  sim.verbose     = 0
  sim.print_detector_info
  sim.set_database(detector_configuration: File.join(script_dir, "../database/detector_configuration.xml"),
                   detector_parameters: detector_parameter_file)
  sim.set_gdml File.join(script_dir, "../database/mass_model_PbCollimator_run11_actual.gdml")
  sim.set_physics(hadron_hp: false, cut_value: 0.001)

  sim.set_primary_generator :PointSourcePrimaryGen, {
    particle: "gamma",
    spectral_distribution: "gaussian",
    energy_mean: energy,
    energy_sigma: 0.0,
    position: vec(32.0+10.0+14.5, 0.0, -12.5),
    direction: vec(-1.0, 0.0, 0.0),
    theta_min: 0.0,
    theta_max: 180.0*Math::PI/180.0,
  }

  sim.run(num)
end

### main ###
random    = 0
script_dir = File.expand_path(__dir__)
num       = ARGV[0].to_i
energy    = ARGV[1].to_f
detector_parameter_file = ARGV[2] || File.join(script_dir, "../database/detector_parameters.xml")
detector_parameter_file = File.expand_path(detector_parameter_file)
random    = (ARGV[3] || random).to_i
output    = ARGV[4] || "simulation_#{num}_#{energy.to_i}keV.root"
p output
run_simulation(num, energy, detector_parameter_file, random, output)
