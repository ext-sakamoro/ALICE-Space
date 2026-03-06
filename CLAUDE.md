# ALICE-Space — Claude Code 設定

## プロジェクト概要

Deep-space communication autonomous trajectory control model-differential protocol

| 項目 | 値 |
|------|-----|
| クレート名 | `alice-space` |
| バージョン | 0.1.0 |
| ライセンス | AGPL-3.0-only |
| リポジトリ | `ext-sakamoro/ALICE-Space` |
| Eco-Systemブリッジ | `bridge_space.rs` + `bridge_space_cross.rs` |
| features | `ffi` |

## 品質スコア

| チェック | 結果 |
|---------|------|
| テスト | 145 passed (144 unit + 1 doc-test) |
| clippy | 0 warnings (`-D warnings`) |
| fmt | clean |
| doc | 0 warnings |
| cdylib | `crate-type = ["rlib", "cdylib"]` |
| FFI | 12関数 (`src/ffi.rs`, feature `ffi`) |
| Unity | `bindings/unity/AliceSpace.cs` |
| UE5 | `bindings/ue5/AliceSpace.h` |

## コーディングルール

メインCLAUDE.md「Git Commit設定」参照。日本語コミット・コメント、署名禁止、作成者 `Moroya Sakamoto`。

## ALICE 品質基準

ALICE-KARIKARI.md「100/100品質基準」参照。clippy基準: `pedantic+nursery`

| 指標 | 値 |
|------|-----|
| clippy (pedantic+nursery) | 0 warnings |
| テスト数 | 145 |
| fmt | clean |

## Eco-System パイプライン

本クレートはALICE-Eco-Systemの以下のパスで使用:
- Path O (Deep-Space→Edge→Analytics/DB/Cache)

## 情報更新ルール

- バージョンアップ時: このCLAUDE.mdのバージョンを更新
- APIの破壊的変更時: ALICE-Eco-Systemブリッジへの影響をメモ
- テスト数/品質の変化時: 品質基準セクションを更新
- 新feature追加時: プロジェクト概要テーブルを更新
