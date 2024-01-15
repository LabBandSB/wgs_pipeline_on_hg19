to run msap you first need to insert sample list into samples.txt

i.e.
# samples.txt
KAZ_WG_001
KAZ_WG_002
KAZ_WG_003

then
python generate_msap_sh.py -j settings.json -s samples.txt

as the result in scripts file will be created new dir named MSAP with script for multi-sample alignment pipeline

