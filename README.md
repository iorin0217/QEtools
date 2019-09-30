# QEtools
QEの計算を便利にするツール
## ディレクトリ
細かい説明はそれぞれ参照
### calc
[pymatgen](http://pymatgen.org/)を使ってcifファイルの情報を読み込み、[ASE](https://wiki.fysik.dtu.dk/ase/)経由で、QEを実行する。
### plot
plotlyを使って以下の項目をプロットする。
- バンド&DOSと、物理量の射影
- フェルミ面と、物理量の射影
- 結晶と、物理量の実空間射影
## ファイル
### auto_qe.py
ジョブを流す
### config.json
パス等の設定
## requirements
- pipenv
## 参考
[pyprocar](https://romerogroup.github.io/pyprocar/)と[ase-espresso](https://ase-espresso.readthedocs.io/en/latest/index.html)
現状とりあえずase-espressoのコードをそのまま使ってるのでGNUライセンス