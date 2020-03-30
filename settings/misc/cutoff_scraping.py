#%%
from selenium import webdriver
from selenium.webdriver.support.ui import Select
import chromedriver_binary
from bs4 import BeautifulSoup
import re
import time
import json
from bokeh.sampledata.periodic_table import elements

browser = webdriver.Chrome()
url = "http://www.pseudo-dojo.org/index.html"
browser.get(url)
type_select = browser.find_element_by_css_selector("#TYP")
type_select.click()
#%%
type_select = Select(type_select)
type_select.select_by_value("nc-sr-04")
#%%
xc = browser.find_element_by_css_selector("#XCF")
xc.click()
#%%
(Select(xc)).select_by_value("pbe")
#%%
acc = browser.find_element_by_css_selector("#ACC")
acc.click()
#%%
(Select(acc)).select_by_value("standard")
data = browser.page_source
html = BeautifulSoup(data)
browser.quit()
#%%
table = html.find_all("div", id=re.compile("_hn"))
cutoff_dict = {}
for div in table[2:]:
    element = div["id"][:-3]
    if not div.text=="na":
        indict = {"cutoff": float(div.text), "dual":4.0, "filename": element+".UPF", "pseudopotential": "ONCVPSP v0.4"}
        cutoff_dict[element] = indict
df = elements.copy()
elements = [i for i in df["symbol"]]
cutoff_dict
#%%
elements
#%%
import os
os.getcwd()
#%%
with open("MIutils/SSSP_efficiency_pseudos/sssp_efficiency.json", "r") as f:
    sssp_dict = json.load(f)
sssp_dict
#%%
with open("MIutils/ONCV_PBE_sr/oncv_pbe_sr.json", "w") as f:
    json.dump(cutoff_dict, f)
with open("MIutils/misc/elements.txt", "w") as f:
    print(*elements,sep="\n",end="\n",file=f)

#%%
no_list_oncv=list(set(elements)-set(cutoff_dict.keys()))
with open("MIutils/ONCV_PBE_sr/no_list.txt", "w") as f:
    print(*no_list_oncv,sep="\n",end="\n",file=f)
no_list_sssp=list(set(elements)-set(sssp_dict.keys()))
with open("MIutils/SSSP_efficiency_pseudos/no_list.txt", "w") as f:
    print(*no_list_oncv,sep="\n",end="\n",file=f)

#%%
