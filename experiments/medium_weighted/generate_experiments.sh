#!/bin/bash
###############################################################################
# generate_experiments.sh
# #
# Generates commands for the different experiments as a file structure in the
# 'generated' folder. The experiments are generated for every possible
# combination of the provided arguments.
# In addition, for every algorithm, the file 'all_commands.txt' is created,
# that holds all commands of the created experiments for the given algorithm.
# #
# Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
###############################################################################

# Import utility methods
source .././utils.sh

# Paths
PATH_TO_HYPERGRAPHS=`cat paths/path_to_hypergraphs.txt`;
PATH_TO_GENERATED_EXPERIMENTS_FOLDER=`cat paths/path_to_generated_experiments_folder.txt`;
PATH_TO_ALGORITHMS=`cat paths/path_to_algorithms.txt`;

# List of algorithms
ALGORITHMS="kernelizer submodular ilp" # Possible algorithms: ilp, kernelizer, submodular, trimmer, kcore_generator

# List of number of label propagation iterations (only for kernelizer and kernelizer_parallel)
LP_NUM_ITERATIONS="0 1";

# List of seeds
SEEDS="0";

# Underlying base solver after kernelization
BASE_SOLVER="submodular"; # Possible solvers: ilp, submodular

# Pruning mode for the kernelizer
PRUNING_MODE="best"; # Possible modes: best, all

# FILE format of the hypergraphs
FILE_FORMAT="HMETIS"; # Possible formats: HMETIS, METIS

# Preset type when reading the hypergraphs
PRESET_TYPE="DETERMINISTIC"; # Possible types: DETERMINISTIC, DEFAULT

# Ordering type for the submodular solver
ORDERING_TYPE="tight" # Possible types: mix_uniform (if parallel), mix_discrete (if parallel), max_adjacency, tight, queyranne

# Ordering mode for the submodular solver
ORDERING_MODE="single" # Possible modes: single, multi

# Number of threads used (only by the parallel algorithms)
# 0 means all available threads
NUM_THREADS=1;

# Memory limit for each instance (in KB)
GB_MEMORY_LIMIT=100000000;

# Set to 1 if you wish to process the hypergraph as unweighted
UNWEIGHTED=0   # set to 1 to enable
UNWEIGHTED_FLAG=""
if [ "$UNWEIGHTED" -eq 1 ]; then
    UNWEIGHTED_FLAG="--unweighted"
fi

# Clear old generated experiments
if [ -d "$PATH_TO_GENERATED_EXPERIMENTS_FOLDER" ]
then
    rm -r ${PATH_TO_GENERATED_EXPERIMENTS_FOLDER};
fi

# Loop over all hypergraphs
for path_to_hypergraph in `cat ${PATH_TO_HYPERGRAPHS}`
do
    hypergraph_name=`basename $path_to_hypergraph | sed 's/\.hgr.*//'`;
    
    # Loop through all seeds
    for seed in $SEEDS
    do
        # Build and create the path where the current experiment should be stored
        path="${PATH_TO_GENERATED_EXPERIMENTS_FOLDER}/jobs/${hypergraph_name}/seed_${seed}";
        
        
        ######################################## trimmer ########################################
        
        if  exists_in_list "${ALGORITHMS}" trimmer
        then
            create_folders_for_algorithm_and_path $path trimmer;
            trimmer_command="${PATH_TO_ALGORITHMS}/heicut_trimmer ${path_to_hypergraph} --seed=${seed} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} --ordering_type=${ORDERING_TYPE} --ordering_mode=${ORDERING_MODE} ${UNWEIGHTED_FLAG}";
            store_command $path trimmer "${trimmer_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## submodular ########################################
        
        if  exists_in_list "${ALGORITHMS}" submodular
        then
            submodular_algorithm_name="submodular_${ORDERING_TYPE}_${ORDERING_MODE}"
            create_folders_for_algorithm_and_path $path $submodular_algorithm_name;
            submodular_command="${PATH_TO_ALGORITHMS}/heicut_submodular ${path_to_hypergraph} --seed=${seed} --file_format=${FILE_FORMAT} --ordering_type=${ORDERING_TYPE} --ordering_mode=${ORDERING_MODE} ${UNWEIGHTED_FLAG}";
            store_command $path $submodular_algorithm_name "${submodular_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## submodular_parallel ########################################
        
        if  exists_in_list "${ALGORITHMS}" submodular_parallel
        then
            submodular_parallel_algorithm_name="submodular_parallel_${ORDERING_TYPE}_${ORDERING_MODE}_T${NUM_THREADS}"
            create_folders_for_algorithm_and_path $path $submodular_parallel_algorithm_name;
            submodular_parallel_command="${PATH_TO_ALGORITHMS}/heicut_submodular_parallel ${path_to_hypergraph} --seed=${seed} --file_format=${FILE_FORMAT} --ordering_type=${ORDERING_TYPE} --ordering_mode=${ORDERING_MODE} --num_threads=${NUM_THREADS} ${UNWEIGHTED_FLAG}";
            store_command $path $submodular_parallel_algorithm_name "${submodular_parallel_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## ilp ########################################
        
        if  exists_in_list "${ALGORITHMS}" ilp
        then
            create_folders_for_algorithm_and_path $path ilp;
            ilp_command="${PATH_TO_ALGORITHMS}/heicut_ilp ${path_to_hypergraph} --seed=${seed} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} ${UNWEIGHTED_FLAG}";
            store_command $path ilp "${ilp_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## ilp_parallel ########################################
        
        if  exists_in_list "${ALGORITHMS}" ilp_parallel
        then
            ilp_parallel_algorithm_name="ilp_parallel_T${NUM_THREADS}"
            create_folders_for_algorithm_and_path $path $ilp_parallel_algorithm_name;
            ilp_parallel_command="${PATH_TO_ALGORITHMS}/heicut_ilp_parallel ${path_to_hypergraph} --seed=${seed} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} --num_threads=${NUM_THREADS} ${UNWEIGHTED_FLAG}";
            store_command $path $ilp_parallel_algorithm_name "${ilp_parallel_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## maxsat ########################################
        
        if  exists_in_list "${ALGORITHMS}" maxsat
        then
            create_folders_for_algorithm_and_path $path maxsat;
            maxsat_command="${PATH_TO_ALGORITHMS}/heicut_maxsat ${path_to_hypergraph} ${path}/maxsat/wcnf.pipe ${path}/maxsat/sol.pipe --seed=${seed} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} ${UNWEIGHTED_FLAG}";
            store_command $path maxsat "${maxsat_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## maxsat_parallel ########################################
        
        if  exists_in_list "${ALGORITHMS}" maxsat_parallel
        then
            maxsat_parallel_algorithm_name="maxsat_parallel_T${NUM_THREADS}"
            create_folders_for_algorithm_and_path $path $maxsat_parallel_algorithm_name;
            maxsat_parallel_command="${PATH_TO_ALGORITHMS}/heicut_maxsat_parallel ${path_to_hypergraph} ${path}/${maxsat_parallel_algorithm_name}/wcnf.pipe ${path}/${maxsat_parallel_algorithm_name}/sol.pipe --seed=${seed} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} --num_threads=${NUM_THREADS} ${UNWEIGHTED_FLAG}";
            store_command $path $maxsat_parallel_algorithm_name "${maxsat_parallel_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## kernelizer ########################################
        
        if  exists_in_list "${ALGORITHMS}" kernelizer
        then
            for lp_num_iterations in $LP_NUM_ITERATIONS
            do
                kernelizer_algorithm_name="kernelizer_IT${lp_num_iterations}"
                create_folders_for_algorithm_and_path $path $kernelizer_algorithm_name;
                kernelizer_command="${PATH_TO_ALGORITHMS}/heicut_kernelizer ${path_to_hypergraph} --seed=${seed} --lp_num_iterations=${lp_num_iterations} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} --base_solver=${BASE_SOLVER} --pruning_mode=${PRUNING_MODE} --ordering_type=${ORDERING_TYPE} --ordering_mode=${ORDERING_MODE} --verbose ${UNWEIGHTED_FLAG}";
                store_command $path $kernelizer_algorithm_name "${kernelizer_command}" $GB_MEMORY_LIMIT;
            done
        fi
        
        ######################################## kernelizer_parallel ########################################
        
        if  exists_in_list "${ALGORITHMS}" kernelizer_parallel
        then
            for lp_num_iterations in $LP_NUM_ITERATIONS
            do
                kernelizer_parallel_algorithm_name="kernelizer_parallel_IT${lp_num_iterations}_T${NUM_THREADS}"
                create_folders_for_algorithm_and_path $path $kernelizer_parallel_algorithm_name;
                kernelizer_parallel_command="${PATH_TO_ALGORITHMS}/heicut_kernelizer_parallel ${path_to_hypergraph} --seed=${seed} --lp_num_iterations=${lp_num_iterations} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} --base_solver=${BASE_SOLVER} --pruning_mode=${PRUNING_MODE} --ordering_type=${ORDERING_TYPE} --ordering_mode=${ORDERING_MODE} --num_threads=${NUM_THREADS} --verbose ${UNWEIGHTED_FLAG}";
                store_command $path $kernelizer_parallel_algorithm_name "${kernelizer_parallel_command}" $GB_MEMORY_LIMIT;
            done
        fi
        
        ######################################## kcore_generator ########################################
        
        if  exists_in_list "${ALGORITHMS}" kcore_generator
        then
            create_folders_for_algorithm_and_path $path kcore_generator;
            kcore_generator_command="${PATH_TO_ALGORITHMS}/heicut_kcore_generator ${path_to_hypergraph} ${path}/kcore_generator/${hypergraph_name} --seed=${seed} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} ${UNWEIGHTED_FLAG}";
            store_command $path kcore_generator "${kcore_generator_command}" $GB_MEMORY_LIMIT;
        fi
        
        ######################################## hypercactus_generator ########################################
        
        if  exists_in_list "${ALGORITHMS}" hypercactus_generator
        then
            create_folders_for_algorithm_and_path $path hypercactus_generator;
            hypercactus_generator_command="${PATH_TO_ALGORITHMS}/heicut_hypercactus_generator ${path_to_hypergraph} ${path}/hypercactus_generator/${hypergraph_name} --seed=${seed} --file_format=${FILE_FORMAT} --preset_type=${PRESET_TYPE} ${UNWEIGHTED_FLAG}";
            store_command $path hypercactus_generator "${hypercactus_generator_command}" $GB_MEMORY_LIMIT;
        fi
    done
done
