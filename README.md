## Converter from HTree to DataTree format

Installation:

1. set HADES ROOT environment:
	source /cvmfs/hades.gsi.de/install/6.12.06/hydra2-4.9w/defall.sh

2. Install DataTree package from https://gitlab.cern.ch/na61-hic/DataTreeQA according to README.md provided with its code

3. Install HTree_to_DT:
	cd HTree_to_DT
	mkdir build
	cd build
	export DATATREE_HOME=/path/to/data/tree/installation
	cmake ..
	make

4. Launch the executable:
	build/HTree_to_DT [input_file_1,input_file_2,...] [output_file] [number_of_events]
	Filelists can be found in the "fileLists" folder
	Do not use the last argument if you want to process all the events.

5. Launch on batch farm:
	./batch/run.sh [path/to/file_list] [number_of_runs]
	Do not use the last argument if you want to process all runs in the list.
	This will launch the convertation in the array of jobs. 
	Output files and logs are stored in the "output" folder.
	
	
