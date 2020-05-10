# 2020/5/4
import json
import tarfile
import zipfile
import shutil
import glob
import requests
from bs4 import BeautifulSoup
import time
import re
from pymatgen.core.periodic_table import get_el_sp
from pymatgen import Element


def predownload():
    # SSSP precision wget
    sssp_cutoff_table_url = "https://www.materialscloud.org/discover/data/discover/sssp/downloads/sssp_precision.json"
    with open("sssp_precision.json", "wb") as f:
        shutil.copyfileobj(requests.get(
            sssp_cutoff_table_url, stream=True).raw, f)
    sssp_pbe_url = "https://www.materialscloud.org/discover/data/discover/sssp/downloads/SSSP_precision_pseudos.tar.gz"
    with open("SSSP_precision_pseudos.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(sssp_pbe_url, stream=True).raw, f)
    with tarfile.open('SSSP_precision_pseudos.tar.gz', 'r') as tf:
        tf.extractall(path='SSSP_precision_pseudos')
    sssp_pbesol_url = "https://www.materialscloud.org/discover/data/discover/sssp/downloads/tarball_PBEsol_1.1/SSSP_PBEsol_pseudos.tar.gz"
    with open("SSSP_PBEsol_pseudos.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(sssp_pbesol_url, stream=True).raw, f)
    with tarfile.open('SSSP_PBEsol_pseudos.tar.gz', 'r') as tf:
        tf.extractall(path='SSSP_PBEsol_pseudos')
    # Dojo wget
    dojo_cutoff_table_url = "https://raw.githubusercontent.com/abinit/pseudo_dojo/master/pseudo_dojo/pseudos/ONCVPSP-PBE-FR-PDv0.4/stringent.djson"
    with open("stringent.djson", "wb") as f:
        shutil.copyfileobj(requests.get(
            dojo_cutoff_table_url, stream=True).raw, f)
    dojo_pbe_sr_url = "http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_stringent_upf.tgz"
    with open("nc-sr-04_pbe_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(dojo_pbe_sr_url, stream=True).raw, f)
    with tarfile.open('nc-sr-04_pbe_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-sr-04_pbe_stringent_upf')
    dojo_pbe_fr_url = "http://www.pseudo-dojo.org/pseudos/nc-fr-04_pbe_stringent_upf.tgz"
    with open("nc-fr-04_pbe_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(dojo_pbe_fr_url, stream=True).raw, f)
    with tarfile.open('nc-fr-04_pbe_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='/nc-fr-04_pbe_stringent_upf')
    dojo_pbesol_sr_url = "http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbesol_stringent_upf.tgz"
    with open("nc-sr-04_pbesol_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(
            dojo_pbesol_sr_url, stream=True).raw, f)
    with tarfile.open('nc-sr-04_pbesol_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-sr-04_pbesol_stringent_upf')
    dojo_pbesol_fr_url = "http://www.pseudo-dojo.org/pseudos/nc-fr-04_pbesol_stringent_upf.tgz"
    with open("nc-fr-04_pbesol_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(
            dojo_pbesol_fr_url, stream=True).raw, f)
    with tarfile.open('nc-fr-04_pbesol_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-fr-04_pbesol_stringent_upf')
    # SG15 wget
    sg15_url = "https://github.com/pipidog/ONCVPSP/archive/master.zip"
    with open("master.zip", "wb") as f:
        shutil.copyfileobj(requests.get(sg15_url, stream=True).raw, f)
    with zipfile.ZipFile('master.zip', 'r') as z:
        z.extractall('master')
    shutil.move("master/ONCVPSP-master/sg15", ".")
    shutil.rmtree("master")
    # GBRV wget
    gbrv_pbe_url = "https://www.physics.rutgers.edu/gbrv/all_pbe_UPF_v1.5.tar.gz"
    with open("all_pbe_UPF_v1.5.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(gbrv_pbe_url, stream=True).raw, f)
    with tarfile.open('all_pbe_UPF_v1.5.tar.gz', 'r') as tf:
        tf.extractall(path='all_pbe_UPF_v1.5')
    gbrv_pbesol_url = "https://www.physics.rutgers.edu/gbrv/all_pbesol_UPF_v1.5.tar.gz"
    with open("all_pbesol_UPF_v1.5.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(gbrv_pbesol_url, stream=True).raw, f)
    with tarfile.open('all_pbesol_UPF_v1.5.tar.gz', 'r') as tf:
        tf.extractall(path='all_pbesol_UPF_v1.5')
    # TODO : pslib wget from GDrive


def pslib_valenve():
    pass


def dojo_valence():
    pass


def sg15_valence():
    pass


def gbrv_valence():
    pass


if __name__ == "__main__":
    # pre definition
    predownload()  # download pseudo potentials
    '''
    pstypes = {'100PAW': "PAW", 'GBRV-1.5': "USPP", 'SG15': "NC", 'GBRV-1.4': "USPP", '031US': "USPP",
               '100US': "USPP", '031PAW': "PAW", 'SG15-1.1': "NC", 'GBRV-1.2': "USPP", 'Dojo': "NC", 'Wentzcovitch': "PAW"}
    Hubbards = {'Sc': {'U': 2.9}, 'Ti': {'U': 4.4}, 'V': {'U': 2.7}, 'Cr': {'U': 3.5}, 'Mn': {'U': 4.0}, 'Fe': {'U': 4.6}, 'Co': {'U': 5.0}, 'Ni': {'U': 5.1}, 'Cu': {'U': 4.0}, 'Zn': {'U': 7.5}, 'Ga': {'U': 3.9}, 'Sn': {'U': 3.5}, 'Nb': {'U': 2.1}, 'Mo': {'U': 2.4}, 'Ta': {'U': 2.0}, 'W': {'U': 2.2}, 'Tc': {'U': 2.7}, 'Ru': {'U': 3.0}, 'Rh': {'U': 3.3}, 'Pd': {'U': 3.6}, 'Ag': {'U': 5.8}, 'Cd': {'U': 2.1}, 'In': {
        'U': 1.9}, 'Re': {'U': 2.4}, 'Os': {'U': 2.6}, 'Ir': {'U': 2.8}, 'Pt': {'U': 3.0}, 'Au': {'U': 4.0}, 'La': {'U': 8.1, 'J': 0.6}, 'Ce': {'U': 7.0, 'J': 0.7}, 'Pr': {'U': 6.5, 'J': 1.0}, 'Nd': {'U': 7.2, 'J': 1.0}, 'Sm': {'U': 7.4, 'J': 1.0}, 'Eu': {'U': 6.4, 'J': 1.0}, 'Gd': {'U': 6.7, 'J': 0.1}, 'Dy': {'U': 5.6, 'J': 0.0}, 'Tm': {'U': 7.0, 'J': 1.0}, 'Yb': {'U': 7.0, 'J': 0.67}, 'Lu': {'U': 4.8, 'J': 0.95}}
    with open("sssp_precision.json", "r") as f:
        ssspp = json.load(f)
    with open("stringent.djson", "r") as f:
        dojo = json.load(f)
    # create elements.json based on sssp
    elements = {}
    for _element in ssspp:
        SOC = "sr"
        print(_element)
        element = get_el_sp(_element)
        pstype = pstypes[ssspp[_element]["pseudopotential"]]
        if (31 <= element.Z <= 34) | (39 <= element.Z <= 52) | (57 <= element.Z <= 84):
            SOC = "fr"
        sssp = {"filename": ssspp[_element]["filename"], "cutoff": ssspp[_element]["cutoff"],
                "dual": ssspp[_element]["dual"], "cite": ssspp[_element]["pseudopotential"], "resource": "SSSP"}
        pseudo = {"pseudopotential": {"PBE": {"sr": {pstype: sssp}}}}
        # collect fr potential
        if SOC == "fr":
            if sssp["cite"] == "Dojo":
                resource = "Dojo"
                cite = "Dojo"
                cutoff = dojo["pseudos_metadata"][_element]["hints"]["high"]["ecut"]
                dual = 4
                ps_file = _element+".rel.oncvpsp.upf"
                dojo_path = glob.glob(
                    f"/home/CMD35/cmd35stud07/QEtools/settings/misc/ONCVPSP-PBE-FR-PDv0.4/{_element}/*.upf")[0]
                shutil.move(
                    dojo_path, f"/home/CMD35/cmd35stud07/QEtools/settings/pseudos/{ps_file}")
            elif "SG15" in sssp["cite"]:
                resource = "ONCVPSP"
                cite = "SG15"
                cutoff = sssp["cutoff"]
                dual = 4
                ps_file = _element+"_ONCV_PBE_fr.upf"
                shutil.move(
                    f"/home/CMD35/cmd35stud07/QEtools/settings/misc/sg15/{ps_file}", f"/home/CMD35/cmd35stud07/QEtools/settings/pseudos/{ps_file}")
            elif (sssp["cite"] == "Wentzcovitch") | ("GBRV" in sssp["cite"]):
                html = requests.get(pslib_url + _element.lower())
                soup = BeautifulSoup(html.text, "html.parser")
                candidates = [i.getText()
                              for i in soup.select("#content-right td")]
                ps_file_indexs = [j for j, k in enumerate(
                    candidates) if f"Pseudopotential type: {pstype} \nFunctional type: PBE\nNon Linear Core Correction\nFull relativistic\n" in k]
                if len(ps_file_indexs) < 1:
                    print(
                        f"I can't find a fr potential of {_element} in PS Library")
                    ps_file = ""
                    cutoff = 0
                    dual = 0
                    cite = ""
                    resource = ""
                else:
                    if len(ps_file_indexs) > 1:
                        # select more semi core
                        ps_file_index = sorted([(len(candidates[j]), j)
                                                for j in ps_file_indexs], reverse=True)[0][1]
                    else:
                        ps_file_index = ps_file_indexs[0]
                    href = "https://www.quantum-espresso.org" + \
                        soup.select(
                            "#content-right td")[ps_file_index].find("a")["href"]
                    ps_text = requests.get(href).text
                    resource = "pslib"
                    cite = "100PAW"
                    # TODO : extract Suggested minimum cutoff for wavefunctions
                    cutoff = sssp["cutoff"]
                    dual = 8
                    ps_file = soup.select(
                        "#content-right td")[ps_file_index].find("a").getText().strip()
                    print(ps_text, sep='\n',
                          file=open(f"/home/CMD35/cmd35stud07/QEtools/settings/pseudos/{ps_file}", 'w'))
            else:
                # sr is taken from pslibrary
                html = requests.get(pslib_url + _element.lower())
                soup = BeautifulSoup(html.text, "html.parser")
                href = "https://www.quantum-espresso.org" + "/upf_files/" + \
                    _element+".rel-" + sssp["filename"][len(_element)+1:]
                if requests.get(href).status_code == 200:
                    ps_text = requests.get(href).text
                    resource = "pslib"
                    cite = sssp["cite"]
                    cutoff = sssp["cutoff"]
                    dual = sssp["dual"]
                    ps_file = _element+".rel-" + \
                        sssp["filename"][len(_element)+1:]
                    print(ps_text, sep='\n',
                          file=open(f"/home/CMD35/cmd35stud07/QEtools/settings/pseudos/{ps_file}", 'w'))
                else:
                    print(
                        f"I can't find a fr potential of {_element} in PS Library")
                    ps_file = ""
                    cutoff = 0
                    dual = 0
                    cite = ""
                    resource = ""
            ps = {"filename": ps_file, "cutoff": cutoff,
                  "dual": dual, "cite": cite, "resource": resource}
            pseudo["pseudopotential"]["PBE"].update({"fr": {pstype: ps}})
        Hubbard = Hubbards[_element] if _element in Hubbards else None
        pseudo.update(
            {"default": {"SOC": SOC, "Hubbard": Hubbard, "pstype": pstype}})
        elements[_element] = pseudo
    with open("/home/CMD35/cmd35stud07/QEtools/settings/elements.json", "w") as f:
        json.dump(elements, f)
    # remove folder and move downloads to misc


    def pslib_download():
        # pslib_ver_dict
        pslib_list_html = requests.get(
            "https://dalcorso.github.io/pslibrary/PP_list.html")
        pslib_list_soup = BeautifulSoup(pslib_list_html.text, "html.parser")
        pslib_list_text = pslib_list_soup.find(class_="highlight").text
        pslib_list_tmp = [i.strip() for i in pslib_list_text.split("\n")]
        pslib_list = [re.sub(':.*_psl.', ':', i)
                    for i in pslib_list_tmp][:-1]  # the last is empty
        pslib_ver_dict = {}
        [pslib_ver_dict.update({i.split(":")[0]: i.split(":")[1]})
        for i in ps_list]
        pslib_url = "https://www.quantum-espresso.org/pseudopotentials/ps-library/"
    '''
