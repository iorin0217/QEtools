# 2020/5/4
import os
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
from pathlib import Path


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
        tf.extractall(path='SSSP_PBE_precision_pseudos')
    sssp_pbesol_url = "https://www.materialscloud.org/discover/data/discover/sssp/downloads/tarball_PBEsol_1.1/SSSP_PBEsol_pseudos.tar.gz"
    with open("SSSP_PBEsol_pseudos.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(sssp_pbesol_url, stream=True).raw, f)
    with tarfile.open('SSSP_PBEsol_pseudos.tar.gz', 'r') as tf:
        tf.extractall(path='SSSP_PBEsol_pseudos')
    shutil.move(
        "SSSP_PBEsol_pseudos/SSSP_PBEsol_pseudos/SSSP_PBEsol_precision_pseudos", "./SSSP_PBEsol_pseudos")
    shutil.rmtree("SSSP_PBEsol_pseudos")
    # Dojo wget
    dojo_cutoff_table_url = "https://raw.githubusercontent.com/abinit/pseudo_dojo/master/pseudo_dojo/pseudos/ONCVPSP-PBE-FR-PDv0.4/stringent.djson"
    print(requests.get(dojo_cutoff_table_url, stream=True).text,
          file=open("stringent.djson", "w"))
    dojo_pbe_sr_url = "http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_stringent_upf.tgz"
    with open("nc-sr-04_pbe_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(dojo_pbe_sr_url, stream=True).raw, f)
    with tarfile.open('nc-sr-04_pbe_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-sr-04_pbe_stringent_upf')
    dojo_pbe_fr_url = "http://www.pseudo-dojo.org/pseudos/nc-fr-04_pbe_stringent_upf.tgz"
    with open("nc-fr-04_pbe_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(dojo_pbe_fr_url, stream=True).raw, f)
    with tarfile.open('nc-fr-04_pbe_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-fr-04_pbe_stringent_upf')
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
    # pslib compile
    # TODO : compile 100 031 script
    pslib_list_html = requests.get(
        "https://dalcorso.github.io/pslibrary/PP_list.html")
    pslib_list_soup = BeautifulSoup(pslib_list_html.text, "html.parser")
    pslib_list_text = pslib_list_soup.find(class_="highlight").text
    pslib_list_tmp = [i.strip() for i in pslib_list_text.split("\n")]
    pslib_list = [re.sub(':.*_psl.', ':', i)
                  for i in pslib_list_tmp][:-1]  # the last is empty
    pslib_ver_dict = {}
    [pslib_ver_dict.update({i.split(":")[0]: "100" if (i.split(":")[1] in ["1.0.0", "1.0.1"]) else "031"})
     for i in pslib_list]
    '''
    # TODO : download CoInMn from QE home page
    pslib_url = "https://www.quantum-espresso.org/pseudopotentials/ps-library/"
    for special in ["co", "in", "mn"]:
        html = requests.get(pslib_url + special.lower())
        soup = BeautifulSoup(html.text, "html.parser")
        candidates = [i.getText()
                        for i in soup.select("#content-right td")]
        ps_file_indexs = [j for j, k in enumerate(
            candidates) if f"Pseudopotential type: {pstype} \nFunctional type: PBE\nNon Linear Core Correction\nFull relativistic\n" in k]
        if len(ps_file_indexs) < 1:
            print(
                f"I can't find a fr potential of {special} in PS Library")
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
            ps_file = soup.select(
                        "#content-right td")[ps_file_index].find("a").getText().strip()
            print(ps_text, sep='\n', file=open(f"./pslib_CoInMn/{ps_file}", 'w'))
    '''


def pslib_read(psfile):
    # UPF ver2
    flag = True
    start = False
    wfc = []
    ecutwfc = None
    with open(psfile, "r") as f:
        while flag:
            line = f.readline()
            if ("cutoff" in line) and ("wavefunctions" in line):
                tmp = (re.search(r"\d+", line)).group()
                if tmp:
                    ecutwfc = int(tmp) if int(tmp) > 0 else None
            if "Valence configuration" in line:
                start = True
            if start:
                if line.split()[0][0] in ["1", "2", "3", "4", "5", "6", "7"]:
                    wfc.append([int(line.split()[0][0]), int(line.split()[2])])
                elif "Generation configuration" in line:
                    flag = False
    return wfc, ecutwfc


def gbrv_read(psfile):
    # UPF ver1
    flag = True
    start = False
    end = False
    wfc = []
    with open(psfile, "r") as f:
        while flag:
            line = f.readline()
            if "<PP_HEADER>" in line:
                start = True
            if start:
                if "Wavefunctions" in line.split()[0]:
                    end = True
            if end:
                if line.split()[0][0] in ["1", "2", "3", "4", "5", "6", "7"]:
                    wfc.append([int(line.split()[0][0]), int(line.split()[1])])
                elif "</PP_HEADER>" in line:
                    flag = False
    return wfc


def oncv_read(psfile, element):
    flag = True
    start = False
    end = False
    tmp = []
    with open(psfile, "r") as f:
        while flag:
            line = f.readline()
            if "PP_INPUTFILE" in line:
                start = True
            if start:
                if element in line.split():
                    nv = int(line.split()[3])
                elif "energy" in line.split():
                    start = False
                    end = True
            if end:
                if line.split()[0] in ["1", "2", "3", "4", "5", "6", "7"]:
                    tmp.append([int(line.split()[0]), int(line.split()[1])])
                elif "PSEUDOPOTENTIAL AND OPTIMIZATION" in line:
                    flag = False
    wfc = tmp[-1 * nv:]
    return wfc


if __name__ == "__main__":
    # pre definition
    # predownload()  # download pseudo potentials
    # os.makedir("pseudos")
    pstypes = {'100PAW': "PAW", 'GBRV-1.5': "US", 'SG15': "ONCV", 'GBRV-1.4': "US", '031US': "US",
               '100US': "US", '031PAW': "PAW", 'SG15-1.1': "ONCV", 'GBRV-1.2': "US", 'Dojo': "ONCV", 'Wentzcovitch': "PAW"}
    Hubbards = {'Sc': {'U': 2.9}, 'Ti': {'U': 4.4}, 'V': {'U': 2.7}, 'Cr': {'U': 3.5}, 'Mn': {'U': 4.0}, 'Fe': {'U': 4.6}, 'Co': {'U': 5.0}, 'Ni': {'U': 5.1}, 'Cu': {'U': 4.0}, 'Zn': {'U': 7.5}, 'Ga': {'U': 3.9}, 'Sn': {'U': 3.5}, 'Nb': {'U': 2.1}, 'Mo': {'U': 2.4}, 'Ta': {'U': 2.0}, 'W': {'U': 2.2}, 'Tc': {'U': 2.7}, 'Ru': {'U': 3.0}, 'Rh': {'U': 3.3}, 'Pd': {'U': 3.6}, 'Ag': {'U': 5.8}, 'Cd': {'U': 2.1}, 'In': {
        'U': 1.9}, 'Re': {'U': 2.4}, 'Os': {'U': 2.6}, 'Ir': {'U': 2.8}, 'Pt': {'U': 3.0}, 'Au': {'U': 4.0}, 'La': {'U': 8.1, 'J': 0.6}, 'Ce': {'U': 7.0, 'J': 0.7}, 'Pr': {'U': 6.5, 'J': 1.0}, 'Nd': {'U': 7.2, 'J': 1.0}, 'Sm': {'U': 7.4, 'J': 1.0}, 'Eu': {'U': 6.4, 'J': 1.0}, 'Gd': {'U': 6.7, 'J': 0.1}, 'Dy': {'U': 5.6, 'J': 0.0}, 'Tm': {'U': 7.0, 'J': 1.0}, 'Yb': {'U': 7.0, 'J': 0.67}, 'Lu': {'U': 4.8, 'J': 0.95}}
    with open("sssp_precision.json", "r") as f:
        ssspp = json.load(f)
    with open("stringent.djson", "r") as f:
        dojo = json.load(f)
    with open("pslib_ver.json", "r") as f:
        pslib_ver_dict = json.load(f)
    # choose pseudo potentials
    elements = {}
    # TODO : Lanthanoides, Co, In, Mn
    for element in without_Lanthanoides:
        SOC = "sr"
        if (31 <= element.Z <= 34) | (39 <= element.Z <= 52) | (57 <= element.Z <= 84):
            SOC = "fr"
        default = {"pstype": pstypes[ssspp[element]["pseudopotential"]],
                   "SOC": SOC, "Hubbard": Hubbards.get(element)}
        pseudopotential = {"PBE": {"sr": {}, "fr": {}},
                           "PBEsol": {"sr": {}, "fr": {}}}
        # PAW
        for func in ["PBE", "PBEsol"]:
            if default["pstype"] == "PAW":
                psfile = [s for s in os.listdir(
                    f"./SSSP_{func}_pseudos") if re.split("[._-]", s)[0].lower() == f"{element}".lower()][0]
                shutil.copy(f"./SSSP_{func}_pseudos/{psfile}", "./pseudos/")
                wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
                ecutwfc = ssspp[element]["cutoff"]
                cite = ssspp[element]["pseudopotential"]
                resource = "SSSP"
            else:
                candidates = [p.name for p in Path(
                    "./pslib" + pslib_ver_dict[element]).glob(f"{element}.{func.lower()}-*-kjpwaw_*")]
                # take longer name which consider more semi-core states
                psfile = max(candidates, key=len)
                shutil.copy(f"./pslib{pslib_ver_dict}/{psfile}", "./pseudos/")
                wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
                cite = pslib_ver_dict[element] + "PAW"
                resource = "PSLibrary" + pslib_ver_dict[element]
            dual = 8.0 if not (element in ["Mn", "Fe", "Co"]) else 12.0
            pseudopotential[func]["sr"]["PAW"] = {
                "filename": f"{psfile}", "cutoff": f"{ecutwfc}", "dual": f"{dual}", "cite": f"{cite}", "resource": f"{resource}"}
            pseudopotential[func]["fr"]["PAW"] = fr_from_sr(
                pseudopotential[func]["sr"]["PAW"])
        # US
        for func in ["PBE", "PBEsol"]:
            # sr
            if default["pstype"] == "US":
                psfile = [s for s in os.listdir(
                    f"./SSSP_{func}_pseudos") if re.split("[._-]", s)[0].lower() == f"{element}".lower()][0]
                shutil.copy(f"./SSSP_{func}_pseudos/{psfile}", "./pseudos/")
                ecutwfc = ssspp[element]["cutoff"]
                cite = ssspp[element]["pseudopotential"]
                resource = "SSSP"
            else:
                candidates = [p.name for p in Path(
                    f"all_{func.lower()}_UPF_v1.5").glob(f"{element.lower()}_*")]
                psfile = candidates[0]  # only one
                shutil.copy(
                    f"all_{func.lower()}_UPF_v1.5/{psfile}", "./pseudos/")
                ecutwfc = ssspp[element]["cutoff"] if default["pstype"] == "PAW" else max(
                    ssspp[element]["cutoff"] / 2, 40)
                cite = "GBRV-1.5"
                resource = "GBRV1.5"
            wfc = gbrv_read(f"./pseudos/{psfile}")
            dual = 8.0 if not (element in ["Mn", "Fe", "Co"]) else 12.0
            pseudopotential[func]["sr"]["US"] = {
                "filename": f"{psfile}", "cutoff": f"{ecutwfc}", "dual": f"{dual}", "cite": f"{cite}", "resource": f"{resource}"}
            # fr
            candidates = [p.name for p in Path(
                "./pslib" + pslib_ver_dict[element]).glob(f"{element}.{func.lower()}-*-rrkjus_*")]
            # take longer name which consider more semi-core states
            psfile = max(candidates, key=len)
            shutil.copy(f"./pslib{pslib_ver_dict}/{psfile}", "./pseudos/")
            wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
            cite = pslib_ver_dict[element] + "US"
            resource = "PSLibrary" + pslib_ver_dict[element]
            pseudopotential[func]["fr"]["US"] = {
                "filename": f"{psfile}", "cutoff": f"{ecutwfc}", "dual": f"{dual}", "cite": f"{cite}", "resource": f"{resource}"}
        # ONCV
        # sr PBE
        dual = 4.0
        if default["pstype"] == "ONCV":
            psfile = ssspp[element]["filename"]
            shutil.copy(f"./SSSP_PBE_pseudos/{psfile}", "./pseudos/")
            wfc = oncv_read(f"./pseudos/{psfile}", element)
            ecutwfc = ssspp[element]["cutoff"]
            cite = ssspp[element]["pseudopotential"]
            resource = "SSSP"
        else:
            candidates = [p.name for p in Path(
                "./sg15").glob(f"{element}_*_sr.upf")]
            psfile = candidates[0]
            shutil.copy(f"./sg15/{psfile}", "./pseudos/")
            wfc = oncv_read(f"./pseudos/{psfile}", element)
            ecutwfc = dojo["pseudos_metadata"][element]["hints"]["high"]["ecut"]
            cite = "SG15"
            resource = "ONCVPSP"
        pseudopotential["PBE"]["sr"]["ONCV"] = {
            "filename": f"{psfile}", "cutoff": f"{ecutwfc}", "dual": f"{dual}", "cite": f"{cite}", "resource": f"{resource}"}
        # fr PBE
        pseudopotential["PBE"]["fr"]["ONCV"] = fr_from_sr(
            pseudopotential["PBE"]["sr"]["ONCV"])
        # sr PBEsol
        if default["pstype"] == "ONCV":
            psfile = [s for s in os.listdir(
                f"./SSSP_PBEsol_pseudos") if re.split("[._-]", s)[0].lower() == f"{element}".lower()][0]
            shutil.copy(f"./SSSP_PBEsol_pseudos/{psfile}", "./pseudos/")
            wfc = oncv_read(psfile, element)
            ecutwfc = ssspp[element]["cutoff"]
            cite = ssspp[element]["pseudopotential"]
            resource = "SSSP"
        else:
            # change name for conflict resolution
            psfile = f"{element}_ONCV_PBEsol-0.4.dojo.upf"
            shutil.copy(
                f"./nc-sr-04_pbesol_stringent_upf/{element}.upf", f"./pseudos/{psfile}")
            wfc = oncv_read(f"./pseudos/{psfile}", element)
            ecutwfc = dojo["pseudos_metadata"][element]["hints"]["high"]["ecut"]
            cite = "Dojo"
            resource = "PseudoDojo0.4"
        pseudopotential["PBEsol"]["sr"]["ONCV"] = {
            "filename": f"{psfile}", "cutoff": f"{ecutwfc}", "dual": f"{dual}", "cite": f"{cite}", "resource": f"{resource}"}
        # fr PBEsol
        psfile = f"{element}_ONCV_PBEsol_fr-0.4.dojo.upf"
        shutil.copy(
            f"./nc-fr-04_pbesol_stringent_upf/{element}.upf", f"./pseudos/{psfile}")
        wfc = oncv_read(f"./pseudos/{psfile}", element)
        ecutwfc = dojo["pseudos_metadata"][element]["hints"]["high"]["ecut"]
        cite = "Dojo"
        resource = "PseudoDojo0.4"
        pseudopotential["PBEsol"]["fr"]["ONCV"] = {
            "filename": f"{psfile}", "cutoff": f"{ecutwfc}", "dual": f"{dual}", "cite": f"{cite}", "resource": f"{resource}"}
        # save to elements
        elements.update(
            {element: {"pseudopotential": pseudopotential, "default": default}})
    with open("./elements.json", "w") as f:
        json.dump(elements, f)
    # remove folder and move downloads to misc
    '''
    shutil.rmtree("SSSP_PBE_pseudos")
    shutil.rmtree("SSSP_PBEsol_pseudos")
    shutil.rmtree("all_pbe_UPF_v1.5")
    shutil.rmtree("all_pbesol_UPF_v1.5")
    shutil.rmtree("nc-fr-04_pbe_stringent_upf")
    shutil.rmtree("nc-fr-04_pbesol_stringent_upf")
    shutil.rmtree("nc-sr-04_pbe_stringent_upf")
    shutil.rmtree("nc-sr-04_pbesol_stringent_upf")
    shutil.rmtree("pslib031")
    shutil.rmtree("pslib100")
    shutil.rmtree("pslib_CoInMn")
    shutil.rmtree("sg15")
    '''