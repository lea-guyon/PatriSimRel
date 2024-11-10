To reproduce the results obtained in the paper, run the script `main.sh` with the following parameters:

## Scenario bilateral descent, ambilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
descent="bilateral"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  \
mf=0.25 \
mm=0.25  
growth_rate=0
sample_size=20  
nbsimu=200  
cores=30  
nameDir="ambilocal_villages"

## Scenario bilateral descent, strict patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
descent="bilateral"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  \
mf=0.5 \
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
mf=0.5 \
mm=0.25  
growth_rate=0
sample_size=20  
nbsimu=100  
cores=40  
nameDir="loose_patrilocal_villages"

## Scenario bilateral descent, strict matrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
descent="bilateral"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  \
mf=0 \
mm=0.5  
growth_rate=0
sample_size=20  
nbsimu=200  
cores=30  
nameDir="strict_matrilocal_villages"

## Scenario bilateral descent, loose matrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
descent="bilateral"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  \
mf=0.25 \
mm=0.5  
growth_rate=0
sample_size=20  
nbsimu=100  
cores=40  
nameDir="loose_matrilocal_villages"

## Scenario strict patrilineal descent, patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150  \
pM=0  \
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0.5 \
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
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0.5 \
mm=0.25 \
sigma=0  
growth_rate=0
sample_size=20  
nbsimu=100   
cores=40  
nameDir="loose_patrilineal_villages"

## Scenario strict matrilineal descent, matrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150  \
pM=0  \
descent="unilineal"  
descent_rule="matrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0 \
mm=0.5 \
sigma=0.1  
growth_rate=0
sample_size=20  
nbsimu=100   
cores=40  
nameDir="strict_matrilineal_villages"

## Scenario loose matrilineal descent, matrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=true  
transmission="full"  
fission_threshold=150  \
pM=0  \
descent="unilineal"  
descent_rule="matrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0.25 \
mm=0.5 \
sigma=0  
growth_rate=0
sample_size=20  
nbsimu=100   
cores=40  
nameDir="loose_matrilineal_villages"

## Scenario strict patrilineal descent to bilateral descent with patrilocal residence
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150  \
pM=0  \
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=1
K=300  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K  
mf=0.5 \
mm=0 \
sigma=0.1  
growth_rate=0
sample_size=20  
nbsimu=100   
cores=40  
nameDir="patrilineal2bilateral"