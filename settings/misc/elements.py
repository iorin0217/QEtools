import json
from pymatgen.core.periodic_table import get_el_sp
from pymatgen import Element
import requests
from bs4 import BeautifulSoup
import re
import time
import shutil

# TODO : SSSP wget
with open("/home/CMD35/cmd35stud07/QEtools/settings/misc/SSSP_precision.json", "r") as f:
    sssps = json.load(f)
# mv SSSP_precision_pseudos/* ../pseudos/

# TODO : Dojo wget
# https://github.com/abinit/pseudo_dojo
# rm other
with open("/home/CMD35/cmd35stud07/QEtools/settings/misc/ONCVPSP-PBE-FR-PDv0.4/stringent.djson", "r") as f:
    dojo = json.load(f)

# TODO : SG15 wget
# https://github.com/pipidog/ONCVPSP/
# rm other

# pslibrary
ps_url = "https://www.quantum-espresso.org/pseudopotentials/ps-library/"

# pre definition
pstypes = {'100PAW': "PAW", 'GBRV-1.5': "USPP", 'SG15': "NC", 'GBRV-1.4': "USPP", '031US': "USPP",
           '100US': "USPP", '031PAW': "PAW", 'SG15-1.1': "NC", 'GBRV-1.2': "USPP", 'Dojo': "NC", 'Wentzcovitch': "PAW"}
Hubbards = {'Sc': {'U': 2.9}, 'Ti': {'U': 4.4}, 'V': {'U': 2.7}, 'Cr': {'U': 3.5}, 'Mn': {'U': 4.0}, 'Fe': {'U': 4.6}, 'Co': {'U': 5.0}, 'Ni': {'U': 5.1}, 'Cu': {'U': 4.0}, 'Zn': {'U': 7.5}, 'Ga': {'U': 3.9}, 'Sn': {'U': 3.5}, 'Nb': {'U': 2.1}, 'Mo': {'U': 2.4}, 'Ta': {'U': 2.0}, 'W': {'U': 2.2}, 'Tc': {'U': 2.7}, 'Ru': {'U': 3.0}, 'Rh': {'U': 3.3}, 'Pd': {'U': 3.6}, 'Ag': {'U': 5.8}, 'Cd': {'U': 2.1}, 'In': {
    'U': 1.9}, 'Re': {'U': 2.4}, 'Os': {'U': 2.6}, 'Ir': {'U': 2.8}, 'Pt': {'U': 3.0}, 'Au': {'U': 4.0}, 'La': {'U': 8.1, 'J': 0.6}, 'Ce': {'U': 7.0, 'J': 0.7}, 'Pr': {'U': 6.5, 'J': 1.0}, 'Nd': {'U': 7.2, 'J': 1.0}, 'Sm': {'U': 7.4, 'J': 1.0}, 'Eu': {'U': 6.4, 'J': 1.0}, 'Gd': {'U': 6.7, 'J': 0.1}, 'Dy': {'U': 5.6, 'J': 0.0}, 'Tm': {'U': 7.0, 'J': 1.0}, 'Yb': {'U': 7.0, 'J': 0.67}, 'Lu': {'U': 4.8, 'J': 0.95}}
elements = {}
SOC = "sr"
# create elements.json based on sssp
for _element in sssps:
    element = get_el_sp(_element)
    pstype = pstypes[sssps[_element]["pseudopotential"]]
    if (31 <= element.Z <= 34) | (39 <= element.Z <= 52) | (57 <= element.Z <= 84):
        SOC = "fr"
    sssp = {"filename": sssps[_element]["filename"], "cutoff": sssps[_element]["cutoff"],
        "dual": sssps[_element]["dual"], "cite": sssps[_element]["pseudopotential"], "resource": "SSSP"}
    pseudo = {"pseudopotential": {"PBE": {"sr": {pstype: sssp}}}}
    # collect fr potential
    if SOC == "fr":
        if sssp["cite"] == "Dojo":
            resource = "Dojo"
            cite = "Dojo"
            cutoff = dojo["pseudos_metadata"][_element]["hints"]["high"]["ecut"]
            dual = 4
            ps_file = _element+".rel.oncvpsp.upf"
            basename = (dojo["pseudos_metadata"]
                        [_element]["basename"]).split(".")[0]
            shutil.move(
                f"/home/CMD35/cmd35stud07/QEtools/settings/misc/ONCVPSP-PBE-FR-PDv0.4/{_element}/{basename}.upf", f"/home/CMD35/cmd35stud07/QEtools/settings/pseudos/{ps_file}")
        elif "SG15" in sssp["cite"]:
            resource = "ONCVPSP"
            cite = "SG15"
            cutoff = sssp["cutoff"]
            dual = 4
            ps_file = _element+"_ONCV_OBE_fr.upf"
            shutil.move(
                f"/home/CMD35/cmd35stud07/QEtools/settings/misc/sg15/{ps_file}", f"/home/CMD35/cmd35stud07/QEtools/settings/pseudos/{ps_file}")
        elif sssp["cite"] == "Wentzcovitch" | "GBRV" in sssp["cite"]:
            html = requests.get(ps_url + _element.lower())
            soup = BeautifulSoup(html.text, "html.parser")
            candidates = [i.getText()
                                    for i in soup.select("#content-right td")]
            ps_file_indexs = [j for j, k in enumerate(
                candidates) if f"Pseudopotential type: {pstype} \nFunctional type: PBE\nNon Linear Core Correction\nFull relativistic\n" in k]
            if len(ps_file_indexs) < 1:
                print(
                    f"I can't find a fr potential of {_element} in PS Library")
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
            html = requests.get(ps_url + _element.lower())
            soup = BeautifulSoup(html.text, "html.parser")
            href = "https://www.quantum-espresso.org" + "/upf_files" + \
                _element+".rel-" + sssp["filename"][len(_element)+1:]
            if requests.get(href).status_code == 200:
                ps_text = requests.get(href).text
                resource = "pslib"
                cite = sssp["cite"]
                cutoff = sssp["cutoff"]
                dual = sssp["dual"]
                ps_file = _element+".rel-" + sssp["filename"][len(_element)+1:]
                print(ps_text, sep='\n',
                    file=open(f"/home/CMD35/cmd35stud07/QEtools/settings/pseudos/{ps_file}", 'w'))
            else:
                print(
                    f"I can't find a fr potential of {_element} in PS Library")
        ps = {"filename": ps_file, "cutoff": cutoff,
            "dual": dual, "cite": cite, "resource": resource}
        pseudo["pseudopotential"]["PBE"].update({"fr": {pstype: ps}})
    Hubbard = Hubbards[_element] if _element in Hubbards else None
    pseudo.update({"default": {"SOC": SOC, "Hubbard": Hubbard, "pstype": pstype})
    elements[_element] = pseudo
with open("/home/CMD35/cmd35stud07/QEtools/settings/misc/SSSP_precision.json", "w") as f:
    json.dump(elements, f)
# TODO : rm
# rm -r SSSP_precision_pseudos
# rm -r ONCVPSP-PBE-FR-PDv0.4
# rm -r sg15
