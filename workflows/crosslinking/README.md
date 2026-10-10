# 架橋エポキシネットワークの構築

溶媒を加える前の架橋モデルを、公開wheelのPAVES APIで生成します。
モノマーのSMILESからpackし、開始点を作成、初期NPT緩和後に
反応→最小化→NVT→NPTを繰り返します。反応後は局所Sageテンプレートで割り当てを更新します。
本体の実装はwheel内のコンパイル済みモジュールにあり、このスクリプトはAPIの実行例です。

リポジトリのルートから環境を準備します。

```bash
conda env create -f workflows/crosslinking/environment.yml
conda activate paves-reactive
pip install dist/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl
cd workflows/crosslinking
python crosslinking.py --project epoxy --total-atoms 5000 \
  --activation-fraction 0.1 --cycles 50 --target-conversion 0.9 \
  --max-reactions-per-cycle 50 --platform CUDA
```

既定モノマーは三官能TMP triglycidyl etherです。`--smiles`か`--smiles-file`で変更できます。
`--total-atoms`は開始点作成時の追加Hも含めた目標で、整数分子数への丸めを記録します。
初期NPTは300 K・1 atm・1 ns、各cycleは最小化、0.1 fs×1000 stepのNVT、
1 fs×1000 stepのNPTです。時間・温度・圧力・反応上限・距離はCLIで変更できます。

反応規則は同じディレクトリの`activation.yaml`と`reactions.yaml`です。
これは中性OHへ活動状態を付けて伝播させる有効モデルで、帯電したカチオン機構そのものを
再現するモデルではありません。反応探索距離から実際の反応速度や障壁は求まりません。
ルールの定義は[YAML_GUIDE.md](YAML_GUIDE.md)を参照してください。

短いCPU動作確認:

```bash
OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python crosslinking.py --project smoke --total-atoms 500 --activation-count 1 \
  --platform CPU --equilibration-steps 2 --cycles 1 \
  --nvt-steps 2 --npt-steps 2 --max-reactions-per-cycle 1 \
  --reaction-cutoff 10 --max-reaction-cutoff 10 --minimization-iterations 100
```

`results/<project>/`に設定、反応履歴、cycleごとの構造、OVITO軌跡と最終構造を保存します。
`final.mol2`には部分電荷と結合・結合次数、`final.structure.json`にはさらに
安定原子ID、活動状態、周期セルを保存します。吸水計算にはJSONを渡します。
周期境界をまたぐ架橋を通常の有限分子用MOL2 importで読み直さないでください。

同名projectがある場合は停止します。明示的な`--force`はそのprojectの結果を削除します。
架橋構築の途中再開には対応していません。保存構造とMD状態は各cycleで出力します。
CPUの短い実行はパイプラインの確認で、目標架橋率や樹脂の平衡を保証しません。

次は[溶媒packと吸水MD](../sorption/README.md)へ進めます。
