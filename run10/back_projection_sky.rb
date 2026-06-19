#!/usr/bin/env ruby

require 'comptonsoft'

Degree = Math::PI / 180.0

class MyApp < ANL::ANLApp
  attr_accessor :inputs, :output

  def setup()
    add_namespace ComptonSoft

    chain :ConstructDetector
    with_parameters(detector_configuration: "database/detector_configuration.xml")
    chain :ReadComptonEventTree
    with_parameters(file_list: inputs)
    chain :BackProjectionSky
    with_parameters(x_min: -180.0,
                    x_max: +180.0,
                    y_min: -90.0,
                    y_max: +90.0,
                    image_center_theta: 90.0,
                    image_center_phi:   0.0,
                    image_yaxis_theta:  0.0,
                    image_yaxis_phi:    0.0,
                    #num_pixel_x: 100,
                    #num_pixel_y: 100,
                    arm: 3.0,
                    num_points: 10000,)
    chain :SaveData
    with_parameters(output: output)
  end
end

data_group_list = ["z0cm", "z4cm", "zm4cm"]

data_group_list.each do |tag|
  app = MyApp.new
  app.inputs = ["products/#{tag}/compton.root"]
  app.output = "products/#{tag}/compton_image_sky.root"
  app.run(:all)
end
