To reproduce the results obtained in the paper, run the script `main.sh` with the following parameters:

## Scenario bilateral descent, strict patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
descent="bilateral"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  \
mf=0.1 \
mm=0  
growth_rate=0
sample_size=20  
nbsimu=200  
cores=30  
nameDir="strict_patrilocal_villages"

## Scenario bilateral descent, loose patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
descent="bilateral"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  \
mf=0.1 \
mm=0.05  
growth_rate=0
sample_size=20  
nbsimu=100  
cores=40  
nameDir="loose_patrilocal_villages"

## Scenario strict patrilineal descent, patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150  \
pM=0  \
friendlyFission="T" \
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0.1 \
mm=0 \
sigma=0.1  
growth_rate=0
sample_size=20  
nbsimu=100   
cores=40  
nameDir="strict_patrilineal_villages"

## Scenario loose patrilineal descent, patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=true  
transmission="full"  
fission_threshold=150  \
pM=0  \
friendlyFission="T" \
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0.1 \
mm=0.05 \
sigma=0  
growth_rate=0
sample_size=20  
nbsimu=100   
cores=40  
nameDir="loose_patrilineal_villages"

## Scenario strict patrilineal descent to bilateral descent with patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150  \
pM=0  \
friendlyFission="T" \
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=1
K=300  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0.1 \
mm=0 \
sigma=0.1  
growth_rate=0
sample_size=20  
nbsimu=100   
cores=40  
nameDir="patrilineal2bilateral"