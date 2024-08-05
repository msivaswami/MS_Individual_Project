import os
import pandas as pd
from io import StringIO

from beamspy import annotation
from beamspy import auxiliary
from beamspy import in_out

from rdkit import Chem

import sys
import tempfile
import numpy as np

import gzip
import sqlite3


def _remove_elements_from_compositions(records, keep, path_nist_database):

    #path_nist_database = os.path.join(params['paths']['path_resources_TPs'], 'libs', 'nist_database.txt')
    nist_database = auxiliary.nist_database_to_pyteomics(path_nist_database)

    elements = [e for e in nist_database if e not in keep]
    for record in records:
        for e in elements:
            if "composition" in record:
                record["composition"].pop(e, None)
            else:
                record.pop(e, None)
    return records


def convert_smiles_inchi_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.MolToInchiKey(mol).split('-')
    return([inchi, inchikey[0], inchikey[1], inchikey[2]])


def fix_msp(msp_pth):

    """
    fix msp files created on Windows where line endings end up causing doubleline spacing
    """

    # replacement strings
    WINDOWS_LINE_ENDING = b'\r\n'
    UNIX_LINE_ENDING = b'\n'

    # relative or absolute file path, e.g.:
    file_path = msp_pth

    with open(file_path, 'rb') as open_file:
        content = open_file.read()
        
    # Windows ➡ Unix
    content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)

    # Unix ➡ Windows
    # content = content.replace(UNIX_LINE_ENDING, WINDOWS_LINE_ENDING)

    with open(file_path, 'wb') as open_file:
        open_file.write(content)
    return


def run_metfrag(path_library, path_msp, path_out, path_config_template,
                path_metfrag, path_beamspy_adducts, assay):

    path_out_orig = path_out
    if len(path_out_orig.split(' ')) > 1:
        path_out_ = tempfile.TemporaryDirectory()
        path_out = path_out_.name

    # Read in msp and extract ms2 info for each grpid
    #extract peak data from msp file
    msms_info = pd.DataFrame()
    fix_msp(path_msp)
    with open(os.path.join(path_msp), "r") as inp_msms:      
        msms_records = inp_msms.read().split("RECORD_TITLE:")
        #msms_records = [record.split("m/z int. rel.int.\n")[1].replace("\n", "\n").strip("\n") for record in
        #                msms_records[1:]]
        def msms_record_extract(msms_record):
            # print(msms_record)
            grpid = msms_record.split("XCMS groupid (grpid):")[1].split('\n')[0].strip(' ')
            # print([msms_record.split("m/z int. rel.int.\n")[1]])
            df = pd.read_csv(StringIO(msms_record.split("m/z int. rel.int.\n")[1]), sep="\t", names=['mz', 'int', 'relint'])
            df.insert(loc=0, column='grpid', value=[grpid] * df.shape[0])
            # print(df)
            return df

        for i in range(1, len(msms_records)):
            msms_info = pd.concat([msms_info, msms_record_extract(msms_records[i])])

        msms_info.loc[:,'grpid'] = msms_info.loc[:,'grpid'].astype('int')

    # Read data to create library for Metfrag
    inp_smiles = pd.read_csv(os.path.join(path_library), delimiter = '\t')
    inp_smiles = inp_smiles.rename(columns = {'prec_mz': 'parentMZ',
                                    'prec_rt': 'parentRT',
                                    'GRPID': 'Identifier',
                                    'INCHI': 'InChI',
                                    'MF': 'MolecularFormula',
                                    'InChIKey.1': 'InChIKey1',
                                    'InChIKey.2': 'InChIKey2',
                                    'InChIKey.3': 'InChIKey3'})
                                    
    inp_smiles.loc[:,'filename'] = ''
    inp_smiles.loc[:,'MonoisotopicMass'] = inp_smiles.loc[:,'parentMZ']

    # print(inp_smiles.columns)

    inp_smiles.loc[:,'InChIKey'] = inp_smiles[["InChIKey1", "InChIKey2", "InChIKey3"]].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    inp_smiles.loc[:,'Name'] = inp_smiles.loc[:,'Identifier']

    #split adduct and formula
    #adducts = {"[M+H]+": 1, "[M+NH4]+": 18, "[M+Na]+": 23,
    #            "[M+K]+": 39, "[M-H]-": -1, "[M+Cl]-": 35, "[M+Acetate]-": 59}

    #update: more flexible approach to creating adducts dictionary - matched input dataframe (duplicates not allowed)
    adducts = pd.read_csv(path_beamspy_adducts, sep = '\t')
    adducts.loc[:,'rounded_mass'] = np.round(adducts.loc[:,'exact_mass'])
    adducts = dict(zip(adducts.label, adducts.rounded_mass))

    def split_adduct_mf(MF, adducts):
        mf = MF.split(" ")[0]
        adduct = MF.split(" ")[1]
        adduct_id = adducts[adduct]
        return(mf, adduct, adduct_id)

    inp_smiles['MolecularFormula'], inp_smiles['adduct'], inp_smiles['adduct_id'] = zip(*inp_smiles.apply(lambda x: split_adduct_mf(x['MolecularFormula'], adducts=adducts), axis = 1))

    # update grpid if processing isomers individually
    # def update_grpid(data):
    #     data.loc[:,'Identifier'] = [str(data.iloc[id]['Identifier']) + '_{}'.format(id+1) for id in range(0,data.shape[0])]
    #     return(data)
    
    ## TODO: check
    # if xeno_sygma_pars == False:
    #     inp_smiles = inp_smiles.groupby('Identifier').apply(lambda x: update_grpid(x))

    def update_grpid_based_on_adduct(data):
        unique_adducts = data['adduct'].unique()
        adduct_map = {adduct: idx for idx, adduct in enumerate(unique_adducts, start=1)}
        data.loc[:, 'Identifier'] = data['Identifier'].astype(str) + '_' + data['adduct'].map(adduct_map).astype(str)
        return data

    inp_smiles = inp_smiles.groupby('Identifier', group_keys=False).apply(update_grpid_based_on_adduct)

    def metfrag_by_grpid(data, msms_info, path_config_template, path_out):
        
        # print("--------DATA-------------")
        # print(data)

        #get groupid as string
        id = data.loc[:,'Identifier'].astype(str).unique()[0]

        # if "3384" not in id:
        #   return

        path_msms_out = os.path.join(str(path_out), "MSMS_matched_{}_{}.msp".format(assay, id))
        # msms_out = msms_info.loc[msms_info.loc[:,'grpid'].astype(str) == id, ['mz', 'int', 'relint']]
        msms_out = msms_info.loc[msms_info.loc[:,'grpid'].astype(str) == id.split("_")[0], ['mz', 'int', 'relint']]

        msms_out.to_csv(path_msms_out, index=False, header=False, sep='\t')

        if data.loc[:,'SMILES'].isnull().values.any() == False and data.loc[:,'MolecularFormula'].isnull().values.any() == False:
            data.loc[:,'filename'] = "MSMS_matched_{}_{}.msp".format(assay, id)

            data.to_csv(os.path.join(path_out, "lib_{}_{}.csv".format(assay, id)), index=False)

            inp_config = open(os.path.join(path_config_template), "r")
            txt_config = inp_config.read()

            path_temp_config = os.path.join(path_out, "config_{}_{}.txt".format(assay, id))

            with open(path_temp_config, "w", newline='\n') as out_config:
                txt_config = txt_config.replace("[PeakListPath]", path_msms_out)
                txt_config = txt_config.replace("[LocalDatabasePath]", path_out + "/lib_{}_{}.csv".format(assay, id))
                txt_config = txt_config.replace("[IonizedPrecursorMass]", data.loc[:,'parentMZ'].astype(str).unique()[0] )
                txt_config = txt_config.replace("[NeutralPrecursorMolecularFormula]", data.loc[:,'MolecularFormula'].astype(str).unique()[0])
                txt_config = txt_config.replace("[SampleName]", "{}_{}".format(assay, id))
                txt_config = txt_config.replace("[PrecursorIonMode]", data.loc[:,'adduct_id'].astype(str).unique()[0].rsplit('.', 1)[0])
                txt_config = txt_config.replace("[ResultsPath]", path_out)
                out_config.write(txt_config)

            if "win" in sys.platform:
                os.system('java -jar {} {}'.format(path_metfrag, path_temp_config))
            else:
                # https://hub.docker.com/r/nfcore/metaboigniter/dockerfile
                # https://github.com/phnmnl/container-metfrag-cli/blob/develop/Dockerfile
                print('metfrag {}'.format(path_temp_config))
                os.system('metfrag {}'.format(path_temp_config))
        else:
            print("Check inputs")
        
        # os.system('java -jar {} {}'.format(path_metfrag, path_temp_config))
        return

    #metfrag based on Identifier
    inp_smiles.groupby('Identifier').apply(lambda x: metfrag_by_grpid(x, 
                                                                      msms_info=msms_info,
                                                                      path_config_template=path_config_template,
                                                                      path_out=path_out))
    
    #move outputs from source to destination directory. /XX preserves existing files in destination
    if "win" in sys.platform:
        os.system('robocopy "{}" "{}" /MIR /XX'.format(path_out, path_out_orig))
    return


def main():

    assay = "HILIC_POS"

    # path_project = "D:/run-mspurity-metfrag/"
    #path_project = "/rds/projects/2015/viantm-01/users/weberrj/run-mspurity-metfrag"
    path_project = r"C:\Users\Muthu\Desktop\Masters Project\FROM_RALF\tools-workflows\run-mspurity-metfrag"
    
    path_input = os.path.join(path_project, "results", "xcms_mspurity")
    path_output = os.path.join(path_project, "results", "metfrag")
    path_output_beamspy = os.path.join(path_project, "results", "beamspy")
    path_resources = os.path.join(path_project, "workflows", "resources")
    path_scripts = os.path.join(path_project, "workflows", "scripts")

    #read in adducts/isotopes database, assign compounds and write out to compound-specific database
    path_peaklist = os.path.join(path_input, "matched_{}.txt".format(assay))
    if "HILIC" in assay:
        path_adducts = os.path.join(path_resources, "beamspy", "adducts_HILIC_beamspy.txt")
    else:
        path_adducts = os.path.join(path_resources, "beamspy", "adducts_LIPIDS_beamspy.txt")

    path_beamspy_db_out = os.path.join(path_output_beamspy, "beamspy_db_{}.sqlite".format(assay))
    path_beamspy_summary = os.path.join(path_output_beamspy, "beamspy_summary_{}.txt".format(assay))
    ppm = 5.0

    #peaklist = in_out.read_peaklist(path_peaklist)
    # Define the mapping dictionary
    mapping = {
        'name': 'name',
        'mz': 'mz',
        'rt': 'rt',
        'intensity': 'intb'  # assuming 'into' is the intensity column
    }

    # Read the peaklist with the correct mapping
    peaklist = in_out.read_peaklist(path_peaklist, mapping=mapping)

    # Drop duplicate rows based on specific columns
    peaklist = peaklist.drop_duplicates(subset=['name', 'mz', 'rt', 'intensity'])
    
    lib_adducts = in_out.read_adducts(path_adducts, "pos")
    path_library_for_metfrag = os.path.join(path_output, 'library_for_metfrag_{}.txt'.format(assay))

    annotation.annotate_compounds(peaklist=peaklist,
                                  lib_adducts=lib_adducts,
                                  ppm=ppm,
                                  db_out=path_beamspy_db_out,
                                  db_name="hmdb_full_v4_0_20200909_v1")

    beamspy_db_name = "hmdb_full_v4_0_20200909_v1.sql.gz"
    path_db_local = os.path.join(path_resources, "beamspy", beamspy_db_name)

    with gzip.GzipFile(path_db_local, mode='rb') as db_dump:
        conn_cpds = sqlite3.connect(":memory:")
        cursor_cpds = conn_cpds.cursor()
        cursor_cpds.executescript(db_dump.read().decode('utf-8'))
        conn_cpds.commit()
        df = pd.read_sql_query('select id, inchi, inchi_key, name as compound_name from hmdb_full_v4_0_20200909_v1', conn_cpds)
        conn_cpds.close()

    df_summary = annotation.summary(in_out.read_peaklist(path_peaklist), path_beamspy_db_out, single_row=False, single_column=False)
    # df_summary.to_csv(path_beamspy_summary, sep="\t", index=False)

    df_summary = pd.merge(df_summary, df, left_on="compound_id", right_on="id")
    df_summary.drop(columns=['id'], inplace=True)
    df_summary[['InChIKey.1', 'InChIKey.2', 'InChIKey.3']] = df_summary['inchi_key'].str.split('-', expand=True)
    df_summary.rename(columns={"inchi": "INCHI"}, inplace=True)

    df_summary.to_csv(path_beamspy_summary, sep="\t", index=False)

    #create Library_Input.csv
    path_peaklist_xlsx = os.path.join(path_input, "matched_{}.xlsx".format(assay))
    msP_matched = pd.read_excel(path_peaklist_xlsx)
    library_input = pd.merge(msP_matched[['grpid', 'name']], df_summary,
                             left_on='name', right_on='name').dropna(subset=['molecular_formula'])

    path_mspurity_grp = os.path.join(path_input, "mspurity_grpeddf_{}.xlsx".format(assay))
    msP_grpeddf = pd.read_excel(path_mspurity_grp)
    msP_grpeddf = msP_grpeddf[['grpid', 'inPurity', 'pid', 'purity_pass_flag', 'precurMtchID', 'precurMtchScan', 'precurMtchRT', 'precurMtchMZ', 'precurMtchPPM']]

    library_input = pd.merge(library_input, msP_grpeddf, on='grpid')
    library_input.loc[:,'MF'] = library_input[['molecular_formula', 'adduct']].apply(lambda x: x.str.cat(sep = ' '), axis = 1)

    #rename columns and subset to keep only those required
    library_input = library_input.rename(columns = {'grpid': 'GRPID', 'precurMtchMZ': 'prec_mz',
                                    'precurMtchRT': 'prec_rt', 'name': 'xcms.id'})[['GRPID', 'prec_mz', 'prec_rt',
                                                                                    'xcms.id', 'MF', 'inPurity',
                                                                                    'INCHI', 'InChIKey.1', 'InChIKey.2', 'InChIKey.3']]
    #                                'compound_name': 'SMILES'})[['GRPID', 'prec_mz', 'prec_rt', 'xcms.id', 'SMILES', 'MF', 'inPurity']]

    library_input["SMILES"] = [Chem.MolToSmiles(Chem.rdinchi.InchiToMol(inchi)[0]) for inchi in library_input["INCHI"]]

    #convert SMILES to INCHI and InChIKey, as per original R script
    library_input.to_csv(path_library_for_metfrag, sep='\t', index=False)

    # DONE ALREADY AGOVE
    # library_input['INCHI'], library_input['InChIKey.1'], library_input['InChIKey.2'], library_input['InChIKey.3'] = zip(*library_input.apply(lambda x: convert_smiles_inchi_inchikey(x['SMILES']), axis = 1))


    #################
    # run metfrag
    #################
    # path_metfrag = os.path.join(path_scripts, "MetFrag2.4.5-CL.jar")
    path_metfrag = os.path.join(path_scripts, "MetFragCommandLine-2.5.0.jar")
    path_msp = os.path.join(path_input, "mspurity_MS2_{}.msp".format(assay))

    path_config = os.path.join(path_resources, "metfrag", "config_template.txt")

    run_metfrag(path_library=path_library_for_metfrag,
                path_msp=path_msp,
                path_out=path_output,
                path_config_template=path_config,
                path_metfrag=path_metfrag,
                path_beamspy_adducts=path_adducts,
                assay=assay)


if __name__ == '__main__':
    sys.exit(main())  # next section explains the use of sys.exit