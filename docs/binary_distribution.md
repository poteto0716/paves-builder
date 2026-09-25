# Binary distribution policy

このリポジトリは polypaves の利用用配布物であり、開発用ソースリポジトリではありません。

## 公開するもの

- `dist/` の CPython/Linux wheel
- Python API の利用例と公開ドキュメント
- 実行に必要な force-field descriptor と明記された第三者データ
- ライセンス、第三者 notice、checksum

## 公開しないもの

- C/C++ の実装、ヘッダー、CMake 入力
- Python 高水準 API の元ソース
- object/static library、ビルドツリー、compile database
- debug symbol、開発ログ、私有パス、ソース管理情報

wheel には `polypaves._native` と `polypaves.api` の2つの strip 済み拡張を収録します。
`polypaves/__init__.py` は公開 API 名を再 export する短い facade だけです。`api.py`、
生成 C、C++ ソースは wheel に含めません。ビルド時のパスはコンパイラの prefix-map
で除去し、ELF の debug section、通常の symbol table、RPATH/RUNPATH を検査します。
動的シンボルは Python がロードに必要とする `PyInit_api` / `PyInit__native` だけを
公開し、C++ エンジンのクラス名・関数名は export しません。

## 限界

公開された実行可能コードは、形式にかかわらず解析の対象になり得ます。難読化や strip
は解析コストを上げ、偶発的な実装情報の露出を減らしますが、完全な防止策ではありません。
そのため技術的対策だけに依存せず、[LICENSE](../LICENSE) で reverse engineering、
decompile、disassemble、非公開実装の復元を禁止しています。

秘密性をさらに優先する配布では、GitHub に wheel 自体を置かず、認証付き artifact
配布またはサーバー側 API とする必要があります。本リポジトリはローカル実行を可能に
するため wheel を公開する設計です。
