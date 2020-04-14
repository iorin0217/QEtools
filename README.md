# QExPy (Quantum Espresso Express with Python)
Quantum Espresso (QE) からの物性計算を誰でも簡単に利用できるようにする

⚠まだ実装途中です
## 特徴
- 原子ごとに用意されたデフォルトの擬ポテンシャル・Hubbard U value・SOCのon/off
- VESTAファイルで定義されたベクトルを初期磁化としてinputできる
- 簡単な物性の計算までをワンストップでWorkflow化
- 計算結果の見やすいDashBoardと論文用のフォーマットの提供
- ASEやAiiDAと違い、フォルダ分けされたlight weightな実装の実現
### for 実験家
- 環境構築がDockerのインストールのみで終わる
- よく計算する物性に関するWorkflowが定義されている
- cifまたはVESTAファイルのみのinputで自動でWorkflowに沿って物性が計算される
### for モデル計算屋
- QEのoutputがhdf5で提供されているので扱いやすい
- ソースコードが単純で改造しやすい
### for 機械学習屋
- DFTに関する知識がなくても機械学習が行える（要注意）
## 機能
hogehoge
## 未実装項目
- On-The-Fly Machine Learning：Optunaとの接続
- dense grid関連：Wannier90との接続
- Phonon関連：ALAMODEとの接続