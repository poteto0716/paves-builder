# Binary validation

公開 wheel は次を検査対象とします。

- wheel に C/C++/Cython 生成 C、ヘッダー、object、static library がないこと
- 高水準API、reactionサブパッケージ、sorptionと`polypaves._native`がstrip済みELF shared objectであること
- debug section、通常の symbol table、RPATH/RUNPATH、私有ビルドパスがないこと
- 各 native extension の公開シンボルが対応する `PyInit_*` だけであること
- 新規仮想環境へ wheel だけでインストールでき、`import polypaves` が成功すること
- `Polymer`、`Copolymer`、`Solvent`、`System` の公開 API が import できること
- OPLS の小規模系を構築し、座標・型・電荷・結合・組成を参照できること
- examples の全 Python ファイルが構文検査を通ること
- PCFF/IFF の収録データが固定 checksum と一致すること

wheelはCPython ABI固有です。0.4.0はCPython 3.11と3.13 / Linux x86_64用です。
OpenFFを使う架橋構築→水pack→固定濃度MD→完了結果のresumeは3.11で検証します。
3.13では依存ライブラリなしでのimportとネイティブ構築を検証します。
他のPython minor versionやOSを検証済みとは扱いません。

strip やパス検査は秘匿性を完全に保証するものではありません。目的は元ソースと容易な
デバッグ情報を配布しないこと、意図したプラットフォームで API が再現可能に動くことです。

公開前の機械監査は次のコマンドで再実行できます。wheel のビルドスクリプトも同じ監査を
最後に自動実行し、いずれかの条件を満たさなければ失敗します。

```bash
scripts/audit_binary_wheel.sh dist/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl
scripts/audit_binary_wheel.sh dist/polypaves-0.4.0-cp313-cp313-linux_x86_64.whl
```

OpenFF環境へ公開wheelをインストールした後のワークフロー検証:

```bash
python -m unittest discover -s tests -p test_public_reactive_workflows.py -v
```

この検証は本体の開発checkoutを使わず、実装モジュールが`.so`から読み込まれることも確認します。
