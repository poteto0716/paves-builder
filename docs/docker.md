# Docker ComposeでPAVESとOVITO Jupyterを使う

公開wheelを使うUbuntu 24.04／Python 3.11の実行環境です。
PAVESの非公開ソースやビルド環境は不要です。OpenFF／Sage、OpenMM、
MDTraj、OVITOのPythonモジュール、JupyterLabを同じ環境に収録します。
配布wheelに合わせ、コンテナは`linux/amd64`です。

## 起動

Docker EngineとDocker Compose v2を用意し、リポジトリのルートで実行します。
Linuxではコンテナのユーザーをファイルの所有者に合わせます。

```bash
export PAVES_UID=$(id -u)
export PAVES_GID=$(id -g)
docker compose up -d --build
docker compose logs jupyter
```

初回は依存パッケージをダウンロードするため時間とディスク容量が必要です。
検証したイメージのサイズは約5.4 GBです。
ログに表示される`http://127.0.0.1:8888/...token=...`をブラウザで開きます。
JupyterLabの開始ページは`workflows/adhesion/ovito_view.ipynb`です。
ログのURLがコンテナ名を指す場合は、ホスト部分を`127.0.0.1`に置き換えます。
`JUPYTER_PORT`を変更した場合は、URLのポートもその値に置き換えます。

設定は環境変数か、[.env.example](../.env.example)をコピーした`.env`に書けます。
`JUPYTER_PORT`でホスト側のポート、`PAVES_CPU_THREADS`でCPUスレッド数を変更します。
`PAVES_JUPYTER_TOKEN`を空にするとランダムなログイントークンを生成します。
固定トークンを指定する場合も、認証は有効なままです。
ホスト側の公開先は既定で`127.0.0.1`です。別のPCから使うときは、例えば
`ssh -L 8888:127.0.0.1:8888 user@server`で転送して開きます。

```bash
docker compose ps
docker compose down
```

## 可視化ノートブック

`workflows/adhesion/ovito_view.ipynb`を開き、カーネルに
**PAVES (Python 3.11)**を選びます。最初のコードセルの`RUN_DIR`を、表示する
計算結果（例：`projects/pmma/runs/09_pull`）に変更して全セルを実行します。
ノートブックの作業ディレクトリは`workflows/adhesion`です。

このノートブックは密着計算の`final.data`、`topology.pdb`、`trajectory.dcd`を使います。
必要なOVITO用dumpは自動生成します。対話表示ではブラウザのWebGLを使い、
マウス操作、frameスライダー、PNG／MP4の出力を利用できます。
GPUなしのコンテナで画像・動画を保存する場合は、rendererを`Tachyon`にすると
CPUで描画できます。MesaのソフトウェアVulkanドライバーも収録しています。
OVITOの画像出力に必要なドライバーは[公式のインストール説明](https://www.ovito.org/docs/current/python/introduction/installation.html)に従っています。

ノートブックのQt操作に合わせて`ipykernel < 7`を固定しています。
コンテナ内だけでカーネルを7以上に更新しないでください。
OVITO 3.16.1は[公式Python wheel](https://www.ovito.org/docs/current/python/introduction/installation.html)から導入します。
最新のOVITO condaパッケージとOpenFFには共有ライブラリのバージョン衝突があるため、
OVITOのQtライブラリはwheelのPySide6依存関係を使います。
conda-forgeの同名OVITOパッケージはPythonモジュールを含まないため使いません。
実行時にはOVITOの一般的なconda警告やMesa／QtのVulkan診断が表示されることがあります。
本構成ではconda側のQt／PySide6を導入せず、ノートブック・PNG・MP4出力まで検証しています。
デスクトップ版OVITO Proの起動はこのComposeの対象外です。

リポジトリ全体を`/workspace`へbind mountするため、ホストにある計算結果をそのまま
読めます。ノートブックの編集、計算結果、画像、動画もホスト側へ保存され、
`docker compose down`で消えません。イメージには既存の`projects/`や`results/`を含めません。

## コンテナ内で計算する

一回だけコマンドを実行する場合は`run --rm`を使います。
起動中のJupyterコンテナで実行する場合は、`exec`に`micromamba run -n base`を付けて
同じ環境を有効にします。

```bash
docker compose run --rm jupyter python docker/smoke_test.py

# 架橋構築の短いCPU動作確認
docker compose run --rm -w /workspace/workflows/crosslinking jupyter \
  python crosslinking.py --project docker_smoke --total-atoms 500 --activation-count 1 \
  --platform CPU --compression-steps 2 --equilibration-steps 2 --cycles 1 --nvt-steps 2 --npt-steps 2 \
  --max-reactions-per-cycle 1 --reaction-cutoff 10 --max-reaction-cutoff 10 \
  --minimization-iterations 100

# 保存電荷を保持した水の配置
docker compose run --rm jupyter python workflows/sorption/sorption.py pack \
  --network workflows/crosslinking/results/docker_smoke/final.structure.json \
  --smiles O --name water --molecules 3 --cell-scale 1.1 \
  --output workflows/sorption/results/docker_smoke/packed

# 起動中のコンテナで吸水MDを実行
docker compose exec jupyter micromamba run -n base \
  python workflows/sorption/sorption.py simulate \
  --packed workflows/sorption/results/docker_smoke/packed \
  --output workflows/sorption/results/docker_smoke/md --platform CPU
```

本計算の条件は[架橋構築](../workflows/crosslinking/README.md)、
[溶媒pack・吸水MD](../workflows/sorption/README.md)、
[密着計算](../workflows/adhesion/README.md)を参照してください。
架橋／吸水のCLIは既定がCUDAなので、CPUでは`--platform CPU`を指定します。
上の短い架橋計算は動作確認用です。指定した水の量は平衡吸水率を意味しません。
架橋／吸水のframe形式は密着計算のDCDノートブックと異なるため、同ノートブックへ
そのまま渡せません。

## NVIDIA GPUを使う

ホストにCUDA 12.9に対応するNVIDIAドライバーとNVIDIA Container Toolkitを用意します。
[ComposeのGPU割り当て](https://docs.docker.com/compose/how-tos/gpu-support/)を追加した設定を使います。

```bash
docker compose -f compose.yaml -f compose.gpu.yaml up -d --build
docker compose -f compose.yaml -f compose.gpu.yaml run --rm jupyter \
  python docker/smoke_test.py --gpu
```

計算コマンドも同じ2つのComposeファイルを指定し、`--platform CUDA`を使います。
CPUとGPUで同じイメージを使います。GPUの予約は追加設定で行うため、
通常の`docker compose up`はGPUなしで起動できます。

## 検証

`docker/smoke_test.py`は配布実装が`.so`から読み込まれること、分子pack、
Sageによる水の割り当て、OpenMM context、OVITOのwidget生成とPNG出力を確認します。
`--gpu`では実際のCUDA contextも作成します。
`--notebook`を付けると小さなDCDを作成し、既存ノートブックの全セル、
frame／wrapの操作、TachyonでのPNG／MP4出力も実行します。

```bash
docker compose run --rm jupyter python docker/smoke_test.py --notebook
```

公開ワークフローの統合テストは次で実行できます。

```bash
docker compose run --rm jupyter \
  python -m unittest discover -s tests -p test_public_reactive_workflows.py -v
```
