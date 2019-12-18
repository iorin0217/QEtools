# QEtools
QEの計算を便利にするツール
## ディレクトリ
細かい説明はそれぞれ参照
### espresso
[pymatgen](http://pymatgen.org/)を使ってcifファイルの情報を読み込み、[ASE](https://wiki.fysik.dtu.dk/ase/)経由で、QEを実行する。
### plot
plotlyを使って以下の項目をプロットする。
- バンド&DOSと、物理量の射影
- フェルミ面と、物理量の射影
- 結晶と、物理量の実空間射影
- 最局在ワニエ関数
## ファイル
### auto_qe.py
ジョブを流す
### config.json
パス等の設定
### elements.json
各原子に対応するデフォルトの
- SOCの有無
- Uの値
- pseudo potentialのカットオフとファイル名
が格納されている
## setup
```bash
pipenv install
wget https://www.materialscloud.org/discover/data/discover/sssp/downloads/SSSP_efficiency_pseudos.tar.gz
```
## 参考
[pyprocar](https://romerogroup.github.io/pyprocar/)と[ase-espresso](https://ase-espresso.readthedocs.io/en/latest/index.html)
現状とりあえずase-espressoのコードをそのまま使ってるのでGNUライセンス