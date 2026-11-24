# first, change directory to the famdb library location
cd /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/share/RepeatMasker/Libraries/famdb/

# download the partitions you require from Dfam 3.9. In the below, change the numbers or range inside the square brackets to choose your subsets.
# e.g. to download partitions 0 to 10: [0-10]; or to download partitions 3,5, and 7: [3,5,7]; [0-16] is ALL PARTITIONS
curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.[0-16].h5.gz'

# decompress Dfam 3.9 paritions
gunzip *.gz

# move up to RepeatMasker main directory
cd /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/share/RepeatMasker/

# save the min_init partition as a backup, just in case!
mv /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/share/RepeatMasker/Libraries/famdb/min_init.0.h5 /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/share/RepeatMasker/Libraries/famdb/min_init.0.h5.bak

# Rerun RepeatMasker configuration
perl ./configure -libdir /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/share/RepeatMasker/Libraries/ -trf_prgm /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/bin/trf -rmblast_dir /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/bin -hmmer_dir /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/bin -abblast_dir /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/bin -crossmatch_dir /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/work/pipeline/homo_sapiens/work/conda/env-e8cd3caef48f718c-9c7bf9a5b0cb0567b5457ce9573c6341/bin -default_search_engine rmblast

