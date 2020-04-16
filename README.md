# ✨QExPy (Quantum Espresso Express with Python)
Quantum Espresso (QE) からの物性計算を誰でも簡単に利用できるようにする

⚠まだ実装途中です
## 特徴
- 原子ごとに用意されたデフォルトの擬ポテンシャル・Hubbard U value・SOCのon/off
- VESTAファイルで定義されたベクトルを初期磁化としてinputできる
- 簡単な物性の計算までをワンストップでWorkflow化
- 計算結果の見やすいDashBoardと論文用のフォーマットの提供
- 環境によらない実行のしやすさ
- ASEやAiiDAに比べlight weightな実装
### for 実験家
- 環境構築がDockerのインストールのみで終わる
- cifまたはVESTAファイルのみのinputで自動でWorkflowに沿って物性が計算される
- pythonが書けなくてもIgorでplotできるoutputフォーマット
### for モデル計算屋
- ワニエ関数が簡単に作成できる
- ソースコードが単純で改造しやすい
### for 機械学習屋
- 自前でハイスループット計算ができる
## セットアップ
## コードの簡単な説明
spread sheet
## 未実装項目
- On-The-Fly Machine Learning：Optunaとの接続
- Tight Binding Model：Wannier90との接続
- Anharmonic Phonon：ALAMODEとの接続