seed_num=100

# Define paths
OSLOM_code_path=/Users/abry4213/Downloads/OSLOM2
data_path=/Users/abry4213/data/OCDA/
RH_data_og=$data_path/HCP_Connectome/RH.txt
tol=0.3
num_iters=100

# Iterate over repeat_num from 1 to 100
# for repeat_num in {1..100}; do
for repeat_num in {1..10}; do

    # Create temporary directory
    repeat_wd=${data_path}/OSLOM_seed_validation/OSLOM_seed_${seed_num}_repeat_${repeat_num}/
    mkdir -p $repeat_wd
    cd $repeat_wd

    # Copy the right hemisphere connectivity data to the temporary directory
    cp $RH_data_og .

    # Define output file
    output_file=${data_path}/OSLOM_seed_validation/OSLOM_seed_${seed_num}_repeat_${repeat_num}.txt

    # Check if output file exists before running OSLOM
    if ! test -f $output_file; then
        # Run OSLOM
        $OSLOM_code_path/oslom_undir -f RH.txt -w -r $num_iters -t $tol -seed $seed_num

        # Move/rename important data to base location
        mv tp $output_file

        # Remove unneeded files
        rm -r RH.txt_oslo_files
    fi

    # Remove the temporary directory
    rm -r ${data_path}/OSLOM_seed_validation/OSLOM_seed_${seed_num}_repeat_${repeat_num}/
done
