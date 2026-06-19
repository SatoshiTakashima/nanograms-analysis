#!/usr/bin/env python3
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def run_simulation_energy(source_name, energy_str, num_photons):
    cmd_list = [
        "./run_simulation.rb",
        str(num_photons),
        energy_str,
    ]

    print(f"Start: {source_name}, {energy_str} keV, {num_photons:.0f} photons")

    result = subprocess.run(
        cmd_list,
        check=True,
    )

    print(f"Done : {source_name}, {energy_str} keV")

    return source_name, energy_str, result.returncode

def run_simulation_nucleus(source_name, za_pair, num_photons, seed):
    cmd_list = [
        "./run_simulation_nucleus.rb",
        str(num_photons),
        str(seed),
        str(za_pair[0]),
        str(za_pair[1])
    ]

    print(f"Start: {source_name}, Z={za_pair[0]}, A={za_pair[1]}, {num_photons:.0f} photons, seed={seed}")

    result = subprocess.run(
        cmd_list,
        check=True,
    )

    print(f"Done: {source_name}, Z={za_pair[0]}, A={za_pair[1]}, seed={seed}")

    return source_name, za_pair, result.returncode


if __name__ == "__main__":
    num_workers = 3

    #sim_config_dict = {"Co60": (27,60), "Na22": (11,22)}
    sim_config_dict     = {"Na22": (11,22)}
    num_photons_per_job = 1e6
    num_loop            = 3
    jobs = []

    for source_name, za_pair in sim_config_dict.items():
        for i in range(num_loop):
            jobs.append((source_name, za_pair, num_photons_per_job, i))
        #for energy_str, num_photons in energy_dict.items():
        #    jobs.append((source_name, energy_str, num_photons))

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for params_job in jobs:
            source_name, za_pair, num_photons_per_job, seed = params_job
            futures.append(executor.submit(run_simulation_nucleus, source_name, za_pair, num_photons_per_job, seed)) 

        for future in as_completed(futures):
            source_name, za_pair, returncode = future.result()
            print(f"Finished: {source_name}, Z={za_pair[0]}, A={za_pair[1]}, returncode={returncode}")

"""
#!/usr/bin/env python3
import subprocess

if __name__=="__main__":
    simConfigDict = {"Co60": {"1173.2": 1e6, "1332.5": 1e6}, "Na22": {"511.0": 2e6,  "1275.4": 1e6}}
    
    for sourceName in simConfigDict.keys():
        for energyStr in simConfigDict[sourceName].keys():
            numPhotons = simConfigDict[sourceName][energyStr]
            cmdList = ["./run_simulation.rb", str(numPhotons), energyStr]
            subprocess.run(cmdList)


"""
#sim_config_dict = {
#    "Co60": {
#        "1173.2": 1e6,
#        "1332.5": 1e6,
#    },
#    "Na22": {
#        "511.0": 2e6,
#        "1275.4": 1e6,
#    },
#}

#jobs = []

#for source_name, energy_dict in sim_config_dict.items():
#    for energy_str, num_photons in energy_dict.items():
#        jobs.append((source_name, energy_str, num_photons))

#with ProcessPoolExecutor(max_workers=num_workers) as executor:
#    futures = [
#        executor.submit(run_simulation, source_name, energy_str, num_photons)
#        for source_name, energy_str, num_photons in jobs
#    ]

#    for future in as_completed(futures):
#        source_name, energy_str, returncode = future.result()
#        print(f"Finished: {source_name}, {energy_str} keV, returncode={returncode}")
