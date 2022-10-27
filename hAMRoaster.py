#!/usr/bin/env python
# coding: utf-8

# python 3.7.7, pandas 1.1.3, numpy 1.19.2
#
# In[1]:
import numpy as np
import pandas as pd
import os
import re
import sys
import argparse
import warnings
# want to avoid print warnings with pandas merges that can be ignored
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version',
                    version='1.0 (initial release)')
parser.add_argument(
    '--fargene', help='full path to fARGene output, if included')
parser.add_argument(
    '--shortbred', help='full path to shortBRED output (tsv), if included')
parser.add_argument(
    '--shortbred_map', help='full path to shortBRED mapping file, if included and not using default')
parser.add_argument(
    '--abx_map', help='full path to Abx:drug class mapping file, if included')
parser.add_argument(
    '--db_files', help='Path to ontology index files, exclude "/" on end of path. Expecting "/aro_categories_index.tsv" in the directory pointed to in this param. ', default="./db_files")
parser.add_argument(
    '--AMR_key', help='full path to key with known AMR phenotypes, REQUIRED', required=True)
parser.add_argument(
    '--name', help='an identifier for this analysis run, REQUIRED', required=True)
parser.add_argument(
    '--ham_out', help='output file from hAMRonization (tsv), REQUIRED', required=True)
parser.add_argument(
    '--groupby_sample', help='Should results from the mock community key be examined per sample (True), or as one whole community (False)', required=False, default='False')

# pull args
args = parser.parse_args()

# Emily's Point To Files for Github
# read in card ontologies as new key
card_key = pd.read_csv(f"{args.db_files}/aro_categories_index.tsv", sep='\t')

# read in drug class key data
if args.abx_map:
    abx_key_name = args.AMR_key
    abx_key = pd.read_csv(abx_key_name)
else:
    abx_key = pd.read_csv(f"{args.db_files}/cleaned_drug_class_key.csv")

# point to known data
if args.AMR_key:
    mock_name = args.AMR_key
    mock = pd.read_csv(mock_name)
    mock['Sample'] = mock['Sample'].str.lower()
# point to observed data, e.g. hAMRonization output
raw_name = args.ham_out
raw_ham = pd.read_csv(raw_name, error_bad_lines=False, sep="\t")
raw_ham['input_file_name'] = raw_ham['input_file_name'].str.lower()

# get ham sum and add fargene and shortbred
#resx = ["res_5x/", "res_50x/", "res_100x/"]
this_run_res = args.name
res_name = args.name
outdir = str(args.name) + "/"
dir_cmd = "mkdir " + outdir
os.system(dir_cmd)

## read in shortbred

if args.shortbred:
    shortbred_path = args.shortbred
    fs = os.listdir(shortbred_path)
    for f in fs:
        sb_path = str(shortbred_path + "/" + f)
        shortbred = pd.read_csv(sb_path, sep="\t")
        # hits only greater than zero
        shortbred = shortbred[shortbred['Hits'] > 0]
        namekey = f.rstrip("_shortbred.tsv")
        shortbred['input_file_name'] = str(namekey)
        shortbred['analysis_software_name'] = "shortbred"
        # merge "family" with gene symbol as closest match
        shortbred['gene_symbol'] = shortbred['Family']

        # give meaning to shortbred family ## to have this if in the for loop is less efficient but it prevents me from having to think very hard when i am very sleep deprived so we shall proceed as is 
        if args.shortbred_map:
            shortmap_name = args.shortbred_map
            shortmap = pd.csv_csv(shortmap_name, sep="\t")
        else:
            shortmap = pd.read_csv(
                f'{args.db_files}/ShortBRED_ABR_Metadata.tab', sep="\t")
    
    
        shortbred = shortbred.merge(shortmap, how="left")
        shortbred['drug_class'] = shortbred['Merged.ID']

        # merge shortbred and rawham results
        raw_ham = raw_ham.append(shortbred, ignore_index=True)

# integrate fargene results
# note there there is a discrepancy in that the results folder has mixed caps and lower case
# so run "find . -depth | xargs -n 1 rename 's/(.*)\/([^\/]*)/$1\/\L$2/' {} \;" in the command line to fix it
beta_lactam_models = ['class_a', 'class_b_1_2',
                      'class_b_3', 'class_c', 'class_d_1', 'class_d_2']
subdirs = ['class_a', 'class_b_1_2', 'class_b_3', 'class_c',
           'class_d_1', 'class_d_2', 'qnr', 'tet_efflux', 'tet_rpg', 'tet_enzyme']

if args.fargene:
    for model in subdirs:
        # for res in resx:
        #  model + "/hmmsearchresults/contigs-" + model + "-hmmsearched.out"
    
        fs = os.listdir( str(args.fargene + model + "/hmmsearchresults/"))
        for f in fs:
            f = str(f)
            fpath = str(args.fargene + model + "/hmmsearchresults/" + f)
            try:
              file = pd.read_csv(fpath, engine='python', sep="\t", header=1, skipfooter=9)

  	          # filter to just what this iteration is
          	  # if res == this_run_res: 
              repi = len(f)
	            #print(res, "\t", model, "\t", len(f))
              if model in beta_lactam_models:
                resistance = "BETA-LACTAM"
              elif "tet" in model:
                resistance = "tetracycline"
              elif model == "qnr":
                resistance = "QUINOLONE"
	            #print("fARGene\t", resistance)

              df = pd.DataFrame(
                      {'analysis_software_name': 'fARGene', 'drug_class': [resistance]})

              newdf = pd.DataFrame(np.repeat(df.values, repi, axis=0))
              newdf.columns = df.columns
	            # print(newdf)
              raw_ham = raw_ham.append(newdf, ignore_index=True)

            except ValueError:
              pass
##############################
    # now fargene and shortbred have been added to the raw_ham table


# This was initially a manual cleaning process, but after pouring through the results at one resolution (20x), here is the closest code version of the manual cleaning process. Basically, we want to default to CARD drug classes for accession number. Then, give all the genes with the same gene name the CARD drug class (overriding the HAM drug class). Fill in NAs with the HAM drug class.

## read in above
card_key['AMR Gene Family'] = card_key['AMR Gene Family'].str.lower()
# have top fill in protein accession NAs with 0.00 so that join works
card_key['Protein Accession'] = card_key['Protein Accession'].fillna(101010)

card = card_key.melt(id_vars=['AMR Gene Family', 'Drug Class', 'Resistance Mechanism'],
                     value_vars=['Protein Accession', 'DNA Accession'],
                     var_name="accession_source",
                     value_name="accession")

card['drug_class_card'] = card['Drug Class']

# merge raw ham and card on accession / reference_accession
salted_ham = raw_ham.merge(
    card, left_on="reference_accession", right_on="accession", how="left")

# if gene_sumbol and drug class card NaN, filter out
smol_gene = salted_ham[['gene_symbol', 'drug_class_card']
                       ].drop_duplicates().dropna(thresh=2)
smol_dict = pd.Series(smol_gene.drug_class_card.values,
                      index=smol_gene.gene_symbol).to_dict()

# make def to apply dict to fill in NAs
def curedThatHam(x):
    var = ''
    if pd.isnull(x['drug_class_card']):  # if there is no card match from accession
        gene = x['gene_symbol']
        try:
            var = smol_dict[gene]
        except KeyError:
            var = x['drug_class_card']

    return var

# apply
salted_ham['drug_class_card'] = salted_ham.apply(lambda x: curedThatHam(
    x), axis=1)  # override OG col because check confirmed it was ok

salted_ham.analysis_software_name.unique()
salted_ham['drug_class_ham'] = salted_ham['drug_class']
cured_ham = pd.DataFrame(salted_ham)

def simplify_drug_class(dat):

    drug_class = ""

    if not pd.isna(dat.drug_class_card):
        drug_class = str(dat['drug_class_card'])
    elif not pd.isna(dat.drug_class_ham):
        drug_class = str(dat['drug_class_ham'])
    elif not pd.isna(dat['Drug Class']):
        drug_class = dat['Drug Class']  # from shortbred
    else:
        drug_class = str("unknown")

    dat['drug_class'] = drug_class

    return dat

cured_ham = cured_ham.apply(lambda x: simplify_drug_class(x), axis=1)
cured_ham['drug_class'] = cured_ham['drug_class'].str.lower()
cured_ham['drug_class'] = cured_ham['drug_class'].apply(
    lambda x: (re.split(r";|, |: | and ", x)))

# has the generic abx names to drug class
abx_key['abx_class'] = (abx_key['drug_class']
                        .apply(lambda x: x if type(x) == str else "")
                        .apply(lambda x: ''.join(e for e in x if e.isalnum())))


abx_melted = abx_key.melt(value_vars=['Generic name', 'Brand names'], id_vars=[
                          'abx_class'], var_name="abx_type", value_name="abx")
# currently this is only adding generic names column
#abx_melted[abx_melted['abx_type']=="Brand names"]

abx_melted['abx'] = abx_melted['abx'].str.lower().str.strip()
# abx_melted[abx_melted['abx']=='tazobactam']
# next combine the drug clas to the antibiotic in the mock community
mock['Abx'] = mock['Abx'].astype(str)
mock['Abx_split'] = mock['Abx'].apply(
    lambda x: str(x).split('-') if "-" in x else x)
mock['Abx'] = (mock['Abx']
               .apply(lambda x: x if type(x) == str else "")
               .apply(lambda x: ''.join(e for e in x if e.isalnum())))

mock = mock.explode('Abx_split')
mock['Abx_split'] = mock['Abx_split'].str.lower()
## add drug class info to mock known list
merged = mock.merge(abx_melted, left_on='Abx_split',
                    right_on='abx', how='left')
merged['abx_class'] = merged['abx_class'].str.lower()

# need to go make data tidy for this to work sucessfully
merged['Abx_split'][merged['abx_class'].isna()].unique()

# Now we have clean data! Now we want to create true +/- and false +/- in `cured_ham`.
# first filter merged so it is only the resistant ones
resistant_mock = merged[merged['classification'] == "resistant"]
len(resistant_mock['abx'].unique())  # number if resustant antibiotuics
len(resistant_mock['abx_class'].unique())  # number drug classes resistant

sus_mock = merged[merged['classification'] == "susceptible"]
boolean_list = ~sus_mock.Abx.isin(resistant_mock['Abx'])
filtered_sus = sus_mock[boolean_list]
# only 6 antibiotics that are susceptible in the entire mock community
filtered_sus.abx.unique()
filtered_sus.abx_class.unique()

# filter filtered sus so that drug classes are unique to sus, not in resistant group
# prior had to filter by antibiotic tested bc only know at drug level, not drug class.
boolean2 = ~filtered_sus.abx_class.isin(resistant_mock['abx_class'])
smol_sus = filtered_sus[boolean2]
smol_sus.abx_class.unique()  

cured_ham = pd.DataFrame(cured_ham.explode('drug_class'))
cured_ham['drug_class'] = (cured_ham['drug_class']
                           .apply(lambda x: x if type(x) == str else "")
                           .apply(lambda x: ''.join(e for e in x if e.isalnum())))

# Now we have the `cured_ham` dataset, which contains all our observations and is cleaned so that every gene observation has a drug class assigned to it (unless the drug class is unknown).
# Now we want to assign those true/false -/+ values.
# give false + / - values
resistant_mock['abx_class'] = resistant_mock['abx_class'].str.lower()
ref_abx = resistant_mock['abx_class'].str.lower()
ref_sus_abx = smol_sus['abx_class'].str.lower()
ref_abx_df = resistant_mock

def get_posneg(row):

    if row['drug_class'] in ref_abx:
        return "true_positive"
    # I would rather explicitly have it search that it is not here
    elif row['drug_class'] in resistant_mock['Abx_split']:
        return "true_positive"
    elif row['drug_class'] in ref_sus_abx:
        return "false_positive"
    else:
        return 'unknown'
    return

cured_ham['True_Positive'] = (cured_ham
                              .apply(lambda x: x['drug_class'] in ref_abx, axis=1))

cured_ham_dc = pd.DataFrame(cured_ham['drug_class'])
false_negatives = ref_abx_df[~(ref_abx_df['abx_class']
                             .isin(cured_ham_dc['drug_class']))]
abx_melted['abx_class'] = abx_melted['abx_class'].str.lower()
false_negatives = cured_ham_dc[~(cured_ham_dc['drug_class']
                                 .isin(ref_abx_df['abx_class']))]
# false_negatives.shape ## precursor to smol_sus; using smol_sus moving forward

# sorry that the variables make less sense as I add things. this is directly related to my sleep quality
# and how much I would like to be done with this project.


# This is where Brooke Talbot did some data cleaning in R. So I exported the above, and read in the results from her cleaning below. SHould probs rewrite what she did in python :,(
# The following code block is originally written by Brooke in R, I have rewritten it in python so that we can use one script to do everything.

# ## R to Python: Drug Class Cleaning

# based on brooke r code
# Summarizing values
# R code
# HAM$drug_class <- str_remove(HAM$drug_class, "antibiotic")
# HAM$drug_class <- str_remove(HAM$drug_class, "resistant")
# HAM$drug_class <- str_remove(HAM$drug_class, "resistance")

# python
cured_ham['drugclass_new'] = cured_ham.drug_class.str.replace('antibiotic', '')
cured_ham['drugclass_new'] = cured_ham.drugclass_new.str.replace(
    'resistant', '')
cured_ham['drugclass_new'] = cured_ham.drugclass_new.str.replace(
    'resistance', '')

# R code
# HAM <- HAM %>% mutate(class_new = case_when(drug_class %in% c("amikacin", "kanamycin", "streptomycin","tobramycin", "kanamycin","spectinomycin", "gentamicin", "aminoglycoside") ~ "Aminoglycoside",
#                                            drug_class %in% c("phenicol", "chloramphenicol") ~ "Phenicol",
#                                            drug_class %in% c("quinolone", "fluoroquinolone", "ciprofloxacinir", "fluoroquinolones") ~ "Quinolones and Fluoroquinolones",
#                                            drug_class %in% c("macrolide", "erythromycin", "mls", "azithromycin", "telythromycin") ~ "Macrolide",
#                                            drug_class %in% c("tetracycline", "glycylcycline") ~ "Tetracycline",
#                                            drug_class %in% c("ampicillin", "methicillin", "penicillin", "amoxicillinclavulanicacid") ~ "Penicillin",
#                                            drug_class %in% c("colistin", "polymyxin", "bacitracin", "bicyclomycin") ~ "Polypeptides",
#                                            drug_class %in% c("cephalosporin", "cefoxatin", "ceftriaxone") ~ "Cephalosporin",
#                                            drug_class %in% c("carbapenem", "penem", "meropenem") ~ "Carbapenem",
#                                            drug_class %in% c("unclassified", "efflux", "acid", "unknown", "multidrug", "multidrugputative", "mutationsonrrnagene16s") ~ "Unclassified",
#                                            drug_class %in% c("linezolid", "oxazolidinone") ~ "Oxazolidinone",
#                                            drug_class %in% c("betalactam", "penam", "betalactamase") ~ "Unspecified Betalactam",
#                                            drug_class %in% c("acridinedye") ~ "Acridine dye",
#                                            drug_class %in% c("antibacterialfreefattyacids") ~ "Free fatty acids",
#                                            drug_class %in% c("benzalkoniumchloride", "quaternaryammonium") ~ "Benzalkonium chlorides",
#                                            drug_class %in% c("peptide") ~ "Unspecified peptide",
#                                            drug_class %in% c("nucleoside") ~ "Unspecified nucleoside",
#                                            drug_class %in% c("fusidicacid") ~ "Fucidic acid",
#                                            drug_class %in% c("sulfonamides", "sulfisoxazole", "sulfonamide") ~ "Sufonamides",
# drug_class %in% c("coppersilver") ~ "copper,silver",
#                                            drug_class %in% c("phenicolquinolone") ~ "phenicol,quinolone", TRUE ~ drug_class))  %>%
#               separate(class_new, into = c("d1", "d2"), sep = "([;,/])")


# python
# make new dict with drug class:drug coding that brooke implemented
cleaning_class = {"Aminoglycoside": ["amikacin", "kanamycin", "streptomycin", "tobramycin", "kanamycin", "spectinomycin", "gentamicin", "aminoglycoside", 'aminocoumarin', 'Aph3ib', 'aph3ib'],
                  "Phenicol": ["phenicol", "chloramphenicol"],
                  "Quinolones and Fluoroquinolones": ["quinolone", "fluoroquinolone", "ciprofloxacinir", "fluoroquinolones"],
                  "Macrolide": ["macrolide", "erythromycin", "mls", "azithromycin", "telythromycin"],
                  "Tetracycline": ["tetracycline", "glycylcycline"],
                  "Penicillin": ["ampicillin", "methicillin", "penicillin", "amoxicillinclavulanicacid"],
                  "Polypeptides": ["colistin", "polymyxin", "bacitracin", "bicyclomycin"],
                  "Cephalosporin": ["cephalosporin", "cefoxatin", "ceftriaxone"],
                  "Carbapenem": ["carbapenem", "penem", "meropenem"],
                  "Unclassified": ["unclassified", "efflux", "acid", "unknown", "multidrug", "multidrugputative",
                                   # clean up mess from split
                                   "mutationsonrrnagene16s", 'warning', 'geneismissingfromnotesfilepleaseinformcurator',
                                   'ant2ia', 'aph6id', 'monobactam', 'shv52a', 'rblatem1'],  # again, mess from split
                  "Oxazolidinone": ["linezolid", "oxazolidinone"],
                  # no need to be "unspecificed" to removed that string
                  "Betalactam": ["betalactam", "penam", "betalactamase"],
                  "Acridine dye": ["acridinedye"],
                  "Free fatty acids": ["antibacterialfreefattyacids"],
                  "Benzalkonium chlorides": ["benzalkoniumchloride", "quaternaryammonium"],
                  "Unspecified peptide": ["peptide", " peptide"],
                  "Unspecified nucleoside": ["nucleoside", " nucleoside"],
                  "Fucidic acid": ["fusidicacid"],
                  "Sufonamides": ["sulfonamides", "sulfisoxazole", "sulfonamide"],
                  #"copper,silver": ["coppersilver", 'copper, silver'],
                  # based on google deep dive, these are different gens of each other?
                  "quinolone": ["phenicolquinolone"],
                  "penicillin": ['amoxicillin/clavulanicacid'],
                  # added based on reviews and critiques of metal categories
                  "metals": ['copper,silver', 'copper', 'silver', 'nickel', 'iron', 'gold', 'zinc', 'cobalt', 'chromium', 'zinc', 'lead', 'mercury', 'arsenic', 'cadmium', 'antimony', 'boron', 'barium', 'tungsten'] # based on reviewer feedback, just shove all the metals into one "metal" category 
                  }
# invert dictionary so the value is the new value, key is old value
new_d = {vi: k for k, v in cleaning_class.items() for vi in v}
# replace key value in column with value from dict (drug class new)
cooked_ham = cured_ham.replace({"drugclass_new": new_d})

# deep down, i am crying, because new_d worked partially but not completely
# i thought I took care of you, new_d, how could you do me like this

# R code: this part did not need to be copied to python; the R version made some duplicate rows from the transofrmation but the python version doesn't do this
# Rearranging rows that needed to be transposed
# Dups1 <- HAM %>% filter(drug_class %in% c("coppersilver","phenicolquinolone")) %>%
#   select(-d2) %>% rename(class_new = d1)
# Dups2 <- HAM %>% filter(drug_class %in% c("coppersilver","phenicolquinolone")) %>%
#   select(-d1) %>% rename(class_new = d2)
# dupout <- HAM %>% filter(is.na(d2)) %>% select(-d2) %>% rename(class_new = d1)
## thanks Brooke Talbot for the R code!! 

# ## Back to the OG Python Script
#cooked_ham = cooked_ham.append(resfinder_results, ignore_index = True)
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].str.lower()

# remove those with Unclassified drug class
to_drop = ['unclassified', 'unknown', 'other', 'multidrug', '', 'target', 'efflux',  # some of these are artifacts of cleaning/split
           'mutationsonproteingene', 'mfsefflux', 'genemodulating', 'mateefflux', 'rndefflux', 'chloramphenicolacetyltransferasecat', 'abcefflux', 'otherarg', 'genemodulatingefflux', 'rrnamethyltransferase']
cooked_ham = cooked_ham[~cooked_ham['drugclass_new'].isin(to_drop)]

# remove unspecified string bewcause it will mess up matching
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].astype(str)
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('unspecified ', ""))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('quinolones and fluoroquinolones', "quinolone"))

# some individual abx made it to drug class col, changing manually to drug class.
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('telithromycin', "macrolide"))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('daptomycin', "lipopeptide"))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('cefoxitin', "cephalosporin"))

cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('classcbetalactamase', 'betalactam'))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('classabetalactamase', 'betalactam'))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('classbbetalactamase', 'betalactam'))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('classdbetalactamase', 'betalactam'))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('betalactamalternatename', 'betalactam'))

# aminoglycosideaminocoumarin
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('aminoglycosideaminocoumarin', 'aminoglycoside'))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('aminoglycosidealternatename', 'aminoglycoside'))
# phosphotransferase would be filtered out as other drug class, so ignore here
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('aminoglycosidephosphotransferase', 'aminoglycoside'))
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(lambda x: x.replace(
    'aminoglycosidenucleotidyltransferase', 'aminoglycoside'))  # same w nucleo...
cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(lambda x: x.replace(
    'aminoglycosideacetyltransferase', 'aminoglycoside'))  # same w acetyl..
# cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(lambda x: x.replace('aminoglycosideaminocoumarin', 'aminoglycoside')) ## same w acetyl..


cooked_ham['drugclass_new'] = cooked_ham['drugclass_new'].map(
    lambda x: x.replace('tetracyclineefflux', 'tetracycline'))

# those that fall into "other" drug class are left as the abx
# we can keep triclosan bc some abx map to triclosan

def assign_true_pos(x):
    target = x['drugclass_new']
    if args.groupby_sample == "False":
        if target in (ref_abx.unique()):
            return "true_positive"
        elif target in (smol_sus['abx_class'].unique()):
            return "false_positive"
        else:
            return "unknown"
    elif args.groupby_sample == "True": 
        ## get relevant variables filterd for sample relevant info 
        try:
            samp = str(x['input_file_name'])
            min_res_mock = ref_abx_df[ref_abx_df['Sample' == samp]]
            ref_abx_min = min_res_mock['abx_class']
            boolean2 = ~sus_mock.abx_class.isin(ref_abx_min)
            smol_sus_min = sus_mock[boolean2]
        except KeyError:
            min_res_mock = ref_abx_df
            ref_abx_min = min_res_mock['abx_class']
            boolean2 = ~sus_mock.abx_class.isin(ref_abx_min)
            smol_sus_min = sus_mock[boolean2]
            
        if target in (ref_abx_min.unique()):
            return "true_positive"
        elif target in (smol_sus_min['abx_class'].unique()):
            return "false_positive"
        else:
            return "unknown"
    
if args.groupby_sample == "False":
    cooked_ham['true_positive'] = cooked_ham.apply(lambda x: assign_true_pos(x), axis = 1)  # true is true pos, false is false neg
    combo_counts = cooked_ham.true_positive.value_counts()
    combo_name = outdir + "combo_counts_" + res_name + ".txt"
    combo_counts.to_csv(combo_name)
    print("\n")
    # first need to do some grouping by tool in cooked_ham
    grouped_ham = cooked_ham.groupby(['analysis_software_name'])
    ham_name = outdir + "cooked_ham_w_true_pos_" + res_name + ".csv"
    cooked_ham.to_csv(ham_name)
else:
    cooked_ham['true_positive'] = cooked_ham.apply(lambda x: assign_true_pos(x), axis = 1)  # true is true pos
    combo_counts = cooked_ham.true_positive.value_counts()
    combo_name = outdir + "combo_counts_" + res_name + ".txt"
    combo_counts.to_csv(combo_name)
    print("\n")
    # first need to do some grouping by tool in cooked_ham
    grouped_ham = cooked_ham.groupby(['analysis_software_name', 'input_file_name'])
    ham_name = outdir + "cooked_ham_w_true_pos_" + res_name + ".csv"
    cooked_ham.to_csv(ham_name) 

# Analysis below
grp_abx_results = grouped_ham['drugclass_new'].value_counts().to_frame()
name_grp_results = outdir + "grouped_by_tool_drug_class" + res_name + ".csv"
grp_abx_results.to_csv(name_grp_results)
# grp_abx_results


pos_count = grouped_ham['true_positive'].value_counts(
).to_frame().unstack(1, fill_value=0)
#pos_count.columns = ['analysis_software_name', 'positive_classification', 'count']

# get false neg dataset
# copied above but def need here
abx_melted['abx_class'] = abx_melted['abx_class'].str.lower()
# drop other from abx_meltered because we excluded unknowns
abx_melted = abx_melted[abx_melted.abx_class != 'other']

# get true abx not in our sample
if args.groupby_sample == "False":
  not_in_mock = abx_melted[~abx_melted['abx_class'].isin(ref_abx)]
  # list of what would be true negative values
  mock_negatives = not_in_mock['abx_class'].unique()
  # use smol_sus for only what is KNOWN as negative
else:
  for samp in resistant_mock['Sample'].unique():
      ref_abx_sub = resistant_mock[resistant_mock['Sample'] == samp]['abx_class']
      not_in_mock = abx_melted[~abx_melted['abx_class'].isin(ref_abx_sub)]
      not_in_mock['Sample'] = str(samp)
      mock_negatives = pd.concat([not_in_mock['abx_class'], not_in_mock['Sample']]).unique()

def fetch_negatives(df, tool):
    df = df[df.analysis_software_name == tool]
    # filter for the tool
    not_in_ham = abx_melted[~abx_melted['abx_class'].isin(df['drugclass_new'])]
    return not_in_ham['abx_class'].unique()
    
    

tool_list = list(cooked_ham['analysis_software_name'].unique())
negs = {}
total_negatives = []

if args.groupby_sample == "False":
    for i in tool_list:
        intermediate = fetch_negatives(cooked_ham, i)
    # print(i, "\n",len(intermediate))
    #  print(intermediate)
    # these appear to be real true negatives - these Abx classes exist, and they are not in any of the hams
        name = str(i)
    # join whats not in ham and not in mock for true neg
        negs[name] = intermediate
        for n in intermediate:
            if n not in total_negatives and n not in cooked_ham['drugclass_new']:
                total_negatives.append(n)
else:
    for i in cooked_ham['input_file_name'].unique():
        subham = cooked_ham[cooked_ham['input_file_name'] == i]
        for i2 in tool_list:
            intermediate = fetch_negatives(subham, i2)
            name = str(i2)
            # join whats not in ham and not in mock for true neg
            try:
                negs[name][i] = intermediate ##key tool, sample name, antibiotics
            except KeyError:
                negs[name] = {}
                negs[name][i] = intermediate
            for n in intermediate:
                if n not in total_negatives and n not in cooked_ham['drugclass_new']:
                    total_negatives.append(n)


all_pos_abx = resistant_mock['abx_class'].unique()
neg_abx = smol_sus.abx_class.unique()
neg_count = pd.DataFrame(columns=['tool', "input_file_name", "false-neg", "true-neg"])

if args.groupby_sample == "False":
    for tool in negs:
        tool_falseneg_count = 0
        tool_trueneg_count = 0
        negatives = negs[tool]
        for n in negatives:
            #print(tool, n, "\n")
            # print([n])
            if n in neg_abx:  # neg_abx for known negative, mock_negatives for negative results including unknown susceptibility
                #print("true_negative: ", tool, n)
                tool_trueneg_count += 1

            elif n in all_pos_abx:
                #print("false negative: ", tool, n)
                tool_falseneg_count += 1

        df2 = {'tool': tool, 'false-neg': tool_falseneg_count,
               'true-neg': tool_trueneg_count}
        neg_count = neg_count.append(df2, ignore_index=True)

        print("Total false negatives from ", tool, ": ", tool_falseneg_count)
        print("Total true negatives from ", tool, ": ", tool_trueneg_count)
        print("\n")
elif args.groupby_sample == "True":
    for tool in negs:
      for samp in cooked_ham['input_file_name'].unique():
          tool_falseneg_count = 0
          tool_trueneg_count = 0
          negatives = negs[tool][samp]
          for n in negatives:
              #print(tool, n, "\n")
              # print([n])
              if n in neg_abx:  # neg_abx for known negative, mock_negatives for negative results including unknown susceptibility
                  #print("true_negative: ", tool, n)
                  tool_trueneg_count += 1

              elif n in all_pos_abx:
                  #print("false negative: ", tool, n)
                  tool_falseneg_count += 1

          df2 = {'tool': tool, "input_file_name" : samp, 'false-neg': tool_falseneg_count,
                 'true-neg': tool_trueneg_count}
          neg_count = neg_count.append(df2, ignore_index=True)
      print("Total false negatives from ", tool, " : ", tool_falseneg_count)
      print("Total true negatives from ", tool, ": ", tool_trueneg_count) 
      print("\n")
print("\n\n")

if args.groupby_sample == "True":
  pos_count_tr = pos_count.stack().unstack(level = 1).reset_index()
  counts = pd.merge(pos_count_tr, neg_count,  how = "outer", right_on = ['tool', 'input_file_name'], left_on = ['analysis_software_name', 'input_file_name'] )
elif args.groupby_sample == "False":
  counts = pd.merge(pos_count, neg_count, right_on="tool", how="outer", left_index = True)

# what re counts if we merge all results
pos_count_total = cooked_ham['true_positive'].value_counts().to_frame().unstack(1, fill_value=0)

# negs - drugs not detected but that exist in the world according to the key
not_in_ham = ref_abx[~ref_abx.isin(cooked_ham['drugclass_new'])]

# not_in_ham
negs_in_ham = cooked_ham['drugclass_new'][~cooked_ham['drugclass_new'].isin(ref_abx)]
cooked_ham[cooked_ham['drugclass_new'].isin(neg_abx)]

# get false/true negs
tot_trueneg_count = 0
tot_falseneg_count = 0

if args.groupby_sample == "False":
  for n in total_negatives:
      #print(tool, n, "\n")
      #  print([n])
      if n in neg_abx:  # neg_abx for known negative, mock_negatives for negative results including unknown susceptibility
          #print("true_negative: ", tool, n)
          if n not in cooked_ham['drugclass_new']:
              tot_trueneg_count += 1
              #print("true neg: ", n)
      elif n in all_pos_abx:
          #print("false negative: ", tool, n)
          if n not in cooked_ham['drugclass_new']:
              tot_falseneg_count += 1
              #print('false_neg:', n)

elif args.groupby_sample == "True":
  for sampname in cooked_ham['input_file_name'].unique():
      samp_cooked_ham = cooked_ham.loc[cooked_ham['input_file_name'] == sampname]
      for n in total_negatives:
        #print(tool, n, "\n")
        #  print([n])
          if n in neg_abx:  # neg_abx for known negative, mock_negatives for negative results including unknown susceptibility
            #print("true_negative: ", tool, n)
              if n not in samp_cooked_ham['drugclass_new']:
                  tot_trueneg_count += 1
                  #print("true neg: ", n)
          elif n in all_pos_abx:
            #print("false negative: ", tool, n)
              if n not in samp_cooked_ham['drugclass_new']:
                  tot_falseneg_count += 1
                  #print('false_neg:', n)     
else:
  print('The input for --grouby_sample is ambiguous and hAMRoaster has quit. It should be either "True" or "False".')
  
# write in a check for the cols we expect
if ('true_positive', 'false_positive') in counts.columns:
    pass
else:
    counts[('true_positive', 'false_positive')] = int(0)

if ('true_positive', 'unknown') in counts.columns:
    pass
else:
    counts[('true_positive', 'unknown')] = int(0)

# # Thanksgiving Ham
## try to add group by arg for sample names


# The following table is what all this code is for.
# sensitivity / specificity analysis
## sensitivity = true_positives / (true_positives + false_negatives)
try:
    counts['sensitivity'] = counts[('true_positive', 'true_positive')] / (counts[('true_positive', 'true_positive')] + counts['false-neg'])
except ZeroDivisionError:
    print('can\'t calculate sensitivity because no values detected. Are you sure your data looks right?')

# precision = true positives / false_positives + true_positives
counts['precision'] = counts[('true_positive', 'true_positive')] / (
    counts['true_positive', 'false_positive'] + counts[('true_positive', 'true_positive')])
# specificity = true_negative / (true_negative + false_positi
try:
    counts['specificity'] = counts['true-neg'] / \
        (counts['true-neg'] + counts[('true_positive', 'false_positive')])
except ZeroDivisionError:
    print("Can't calculate specificity because there are no observed true negatives or false positives.")
    counts['specificity'] = int(0)
## accuracy = (true_positive + true_negative) / (true_positive + false_positive + true_negative + false_negative)
counts['accuracy'] = (counts[('true_positive', 'true_positive')] + counts['true-neg']) / (counts[('true_positive',
                                                                                                  'true_positive')] + counts['true_positive', 'false_positive'] + counts['true-neg'] + counts['false-neg'])
## percent unclassified, now not a proportion in response to reviewer feedback 
counts['percent_unclassified'] = (counts[('true_positive', 'unknown')] / (counts[(
    'true_positive', 'true_positive')] + counts[('true_positive', 'false_positive')] + counts[('true_positive', 'unknown')])) * 100


#print("Thanksgiving Ham, ", this_run_res, ": ")
name = outdir + "thanksgiving_ham_final_" + res_name + ".csv"
counts.to_csv(name)
# print(counts) ## print out if interactive; does nothing if command line
print("\nYou have new results! Check it out in ", name, "\nHint: Your main results will be thanksgiving_ham_final_<NAME>.")
#####################################################
# # Canned Ham  # # # # # # # #
# condensing the results
# note that these aren't particularly informative. it's here in case reviewers really want it.

cooked_ham = cooked_ham.assign(condensed_gene=cooked_ham.groupby(
    'analysis_software_name')['input_gene_start'].shift(-1))


def condense_results(df):
    start_col = df['input_gene_start']
    stop_col = df['input_gene_stop']
    next_row = df['condensed_gene']
    keep_val = ""
    condense_val = ""
    unconclusive_val = ""

    if next_row < stop_col:  # if the start of the next gene overlaps with the end of THIS gene
        condense_val = "condense"
    elif stop_col < next_row:
        keep_val = "keep"
    else:
        unconclusive_val = "unconclusive"

    message = str(condense_val + keep_val + unconclusive_val)
    return (message)


cooked_ham['condense_action'] = cooked_ham.apply(lambda x: condense_results(
    x), axis=1)  # does this need to be grouped by analysis software name ?
# cooked_ham['condense_action'].describe()
# looks like there are 1811 instances where AMR genes overlap. Lets remove these and see what happens.

canned_ham = cooked_ham[cooked_ham['condense_action'] != "condense"]
# i am very proud of myself for this name
# its the little things in life that bring me smiles

# positives
grouped_can = canned_ham.groupby(['analysis_software_name'])
pos_count2 = grouped_can['true_positive'].value_counts(
).to_frame().unstack(1, fill_value=0)

# add negatives
negs2 = {}
for i in tool_list:
    intermediate = fetch_negatives(canned_ham, i)
    #print(i, "\n",len(intermediate))
    # these appear to be real true negatives - these Abx classes exist, and they are not in any of the hams
    name = str(i)
    # join whats not in ham and not in mock for true neg
    negs2[name] = intermediate


# negative counts
neg_count2 = pd.DataFrame(columns=['tool', "false-neg", "true-neg"])
for tool in negs2:
    tool_falseneg_count = 0
    tool_trueneg_count = 0
    negatives = negs2[tool]
    for n in negatives:
       # print(tool, n, "\n")
        # print([n])
        if n in neg_abx:
            #print("true_negative: ", tool, n)
            tool_trueneg_count += 1

        elif n in neg_abx:  # only things that are known to be negative
            #print("false negative: ", tool, n)
            tool_falseneg_count += 1

    df3 = {'tool': tool, 'false-neg': tool_falseneg_count,
           'true-neg': tool_trueneg_count}
    neg_count2 = neg_count2.append(df3, ignore_index=True)

    #print("False negatives from ", tool, ": ", tool_falseneg_count)
    #print("Truenegatives from ", tool, ": ", tool_trueneg_count)

# merge this all together
counts2 = pd.merge(pos_count2, neg_count2, right_on="tool",
                   how="outer", left_index=True, right_index=False)
# counts2 ## this spits out an error, but we can ignore it
# write in a check for the cols we expect
if ('true_positive', 'false_positive') in counts2.columns:
    pass
else:
    counts2[('true_positive', 'false_positive')] = int(0)

if ('true_positive', 'unknown') in counts2.columns:
    pass
else:
    counts2[('true_positive', 'unknown')] = int(0)

# sensitivity / specificity analysis
## sensitivity = true_positives / (true_positives + false_negatives)
try:
    counts2['sensitivity'] = counts2[('true_positive', 'true_positive')] / (
        counts2[('true_positive', 'true_positive')] + counts2['false-neg'])
except ZeroDivisionError:
    counts2['sensitivity'] = int(0)
# precision
counts2['precision'] = counts2[('true_positive', 'true_positive')] / (
    counts2['true_positive', 'false_positive'] + counts2[('true_positive', 'true_positive')])
## specificity = true_negative / (true_negative + false_positive)
try:
    counts2['specificity'] = counts2['true-neg'] / \
        (counts2['true-neg'] + counts2[('true_positive', 'false_positive')])
except ZeroDivisionError:
    print("Can't calculate specificity because you have no observed true negatives or false positives. Is this expected?")
    counts2['specificity'] = int(0)
    # if we do true neg (aka if it is not zero:)
## accuracy = (true_positive + true_negative) / (true_positive + false_positive + true_negative)
try:
    counts2['accuracy'] = (counts2[('true_positive', 'true_positive')] + counts2['true-neg']) / (counts2[('true_positive',
                                                                                                          'true_positive')] + counts2['true_positive', 'false_positive'] + counts2['true-neg'] + counts2['false-neg'])
except ZeroDivisionError:
    counts2['accuracy'] = 0
# recall
# true pos / (true pos  + false neg)
#counts2['recall'] = counts2[('true_positive','true_positive')] / (counts2[('true_positive', 'true_positive')] + counts2['false-neg'])
#counts2['recall'] = pd.to_numeric(counts2['recall'])
# F1
## 2 * (precision * recall) / (precision + recall)
#counts2['F1'] = 2 * ( (counts2['precision'] * counts2['recall']) / (counts2['precision'] + counts2['recall']) )
# except ZeroDivisionError:
#    counts['F1'] = 0
counts2['percent_unclassified'] = counts2[('true_positive', 'unknown')] / (counts2[(
    'true_positive', 'true_positive')] + counts2[('true_positive', 'false_positive')] + counts2[('true_positive', 'unknown')])

#print("Canned Ham Thanksgiving Table, ", this_run_res, " :")
# no need to print this
name2 = outdir + "canned_ham_" + res_name + ".csv"
counts2.to_csv(name2)
counts2  # prints if interactive, ignored if not

# combine results of all tools
# count specificity
tot_false_pos = sum(counts[('true_positive', 'false_positive')])
tot_true_pos = sum(counts[('true_positive', 'true_positive')])
tot_false_neg = min(counts['false-neg'])
tot_true_neg = min(counts['true-neg'])
tot_unknown = sum(counts['true_positive', 'unknown'])

#################################################################
# combo stats ###################################################
#################################################################
## note that this is not particularly informative so it's here, #
## but it does not do anything for the user. This will just put #
## some output in the interactive STOUT #########################
#################################################################

## sensitivity = true_positives / (true_positives + false_negatives)
sensitivity = tot_true_pos / (tot_true_pos + tot_false_neg)

# precision
precision = tot_true_pos / (tot_false_pos + tot_true_pos)

## specificity = true_negative / (true_negative + false_positive)
try:
    specificity = tot_true_neg / (tot_true_neg + tot_false_pos)
except ZeroDivisionError:
    specificity = "0"
## accuracy = (true_positive + true_negative) / (true_positive + false_positive + true_negative)
accuracy = (tot_true_pos + tot_true_neg) / (tot_true_pos +
                                            tot_false_pos + tot_true_neg + tot_false_neg)

# recall
# true pos / (true pos  + false neg)
#recall = tot_true_pos / (tot_true_pos + tot_false_neg)

# F1
## 2 * (precision * recall) / (precision + recall)
#
#F1 = 2 * ( precision * recall / (precision + recall) )
# except ZeroDivisionError:
#    counts['F1'] = 0
percent_unclassified = tot_unknown / \
    (tot_true_pos + tot_false_pos + tot_unknown)

# not printing this because I think it will confuse users. 
#print("combo stats (compiled outputs of all tools): ", "\n", "sensitivity: ", sensitivity, "\n specificity",
#      specificity, "\n precision", precision, "\n accuracy", accuracy, "\n")  # recall", recall)
#print(" percent unknown: ", percent_unclassified, "\n")#\

#print("all combined basic counts\ntot_false_pos", tot_false_pos, "\ntot_true_pos", tot_true_pos,
#      "\ntot_false_neg", tot_false_neg, "\ntot_true_neg", tot_false_neg, "\ntot_unknown", tot_unknown)


print("\n\n\nRun Complete. Thank you for using hAMRoaster.")
