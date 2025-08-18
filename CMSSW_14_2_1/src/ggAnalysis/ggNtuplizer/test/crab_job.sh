SAMPLES_FILE=${1:="data_2023_C.txt"}

#input parameters
config_file=run.py
publication=False
site=T3_CH_CERNBOX
DBS=global


YEAR=${2:-"2023"}
ERA=${3:-"C"}
IsDATA=${4:-"true"}


identifier=${YEAR}
if [[ "$IsDATA" == "true" ]]; then
    identifier+="_DATA_"${ERA}
else
    identifier+="_MC"
fi

IsDATA_crab=0
if [[ "$IsDATA" == "true" ]]; then
    IsDATA_crab=1
fi

fil_list=crab_submit_${identifier}
mon_list=crab_monitor_${identifier}

# Clear previous submission/monitoring scripts
> ${fil_list}.sh
> ${mon_list}.sh

if [[ ! -f "$SAMPLES_FILE" ]]; then
  echo "Error: Samples file '$SAMPLES_FILE' not found."
  exit 1
fi

# Read sample names from the input file
mapfile -t sample_data < "$SAMPLES_FILE"

counter=1

for sample in "${sample_data[@]}"; do

    sample=$(echo "$sample" | tr -d '\r')

    echo "Processing: ${sample}"

    primary_dataset=$(echo "$sample" | cut -d'/' -f2)      # Should give: EGamma0
    processed_dataset=$(echo "$sample" | cut -d'/' -f3)    # Should give: Run2023C-22Sep2023_v1-v1
    
    echo "  Primary dataset: '$primary_dataset'"
    echo "  Processed dataset: '$processed_dataset'"
    
    if [[ "$IsDATA" == "true" ]]; then
        # For data: use primary dataset + processed info + counter
        # Clean the processed dataset name (replace - with _)
        processed_clean=$(echo "$processed_dataset" | sed 's/-/_/g')
        label="${primary_dataset}_${processed_clean}_${counter}"
    else
        # For MC: use dataset name + counter
        label="${primary_dataset}_${counter}"
    fi
    
    # Additional cleaning - remove any remaining problematic characters
    label=$(echo "$label" | sed 's/+/plus/g' | cut -c1-100)
    
    echo "  Final label: '$label'"
    
    dataset=$(echo "$sample" | tr -d '\r')  # Remove carriage return if present
    
    # Generate crab config
    ./write_crab_config.sh $label $identifier $config_file $dataset $publication $site $DBS $YEAR $IsDATA_crab $ERA
    echo "  Crab config crabfile_${YEAR}_${label}.py created."
    echo ""
    
    # Add to submission and monitoring scripts
    echo "crab submit -c crabfile_${YEAR}_${label}.py" >> ${fil_list}.sh
    echo "crab status -d crab_${identifier}/crab_crab_${label}/" >> ${mon_list}.sh
    
    ((counter++))
done

echo "All crab configs created successfully."

chmod +x ${fil_list}.sh
chmod +x ${mon_list}.sh

echo "Submission script: ${fil_list}.sh"
echo "Monitoring script: ${mon_list}.sh"
