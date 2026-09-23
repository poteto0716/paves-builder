# Binary validation

公開 wheel は次を検査対象とします。

- wheel に C/C++/Cython 生成 C、ヘッダー、object、static library がないこと
- `polyse.api` と `polyse._native` が strip 済み ELF shared object であること
- debug section、通常の symbol table、RPATH/RUNPATH、私有ビルドパスがないこと
- 各 native extension の公開シンボルが対応する `PyInit_*` だけであること
- 新規仮想環境へ wheel だけでインストールでき、`import polyse` が成功すること
- `Polymer`、`Copolymer`、`Solvent`、`System` の公開 API が import できること
- OPLS の小規模系を構築し、座標・型・電荷・結合・組成を参照できること
- examples の全 Python ファイルが構文検査を通ること
- PCFF/IFF の収録データが固定 checksum と一致すること

wheel は CPython ABI 固有です。現在の artifact は CPython 3.13 / Linux x86_64 用で、
他の Python minor version や OS を検証済みとは扱いません。

strip やパス検査は秘匿性を完全に保証するものではありません。目的は元ソースと容易な
デバッグ情報を配布しないこと、意図したプラットフォームで API が再現可能に動くことです。

公開前の機械監査は次のコマンドで再実行できます。wheel のビルドスクリプトも同じ監査を
最後に自動実行し、いずれかの条件を満たさなければ失敗します。

```bash
scripts/audit_binary_wheel.sh dist/polyse-0.2.0-cp313-cp313-linux_x86_64.whl
```
