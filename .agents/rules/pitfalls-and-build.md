# 実装上の落とし穴とビルド手順（coresynth）

## 環境・R・C++ の落とし穴
- **`arma::solve_opts::fast` は使用禁止** — 本環境（Windows + Rtools45）でクラッシュする。
- Armadillo の `vec` は R 側で N×1 行列になる → R 層で `as.numeric()` により正規化すること。
- `paste0("Donor ", integer(0))` は zero-length が `""` 扱いになり、長さ1の `"Donor "` を返してしまう → zero-length を保持する `sprintf()` を使う。
- **`R CMD INSTALL` 失敗時は旧版にサイレントロールバックされる** — 失敗に気づかず古い DLL でテストしてしまう罠がある。ビルド後は必ずインストールが成功したか確認すること。
- `gh` CLI はインストール済み・認証済み（アカウント `yo5uke`）。PR 作成等は `gh` コマンドで行ってよい。
- **plot のタイトル・サブタイトルに非 ASCII 文字を入れない**（ギリシャ文字・Unicode マイナス記号 `−` U+2212 等）。UTF-8 でないロケール（GitHub Actions macOS ランナー等）で grid のテキストメトリクス計算が `mbcsToSbcs` 変換に失敗し、`R CMD check --run-donttest` がクラッシュする。マイナス記号は ASCII ハイフン `-`、ギリシャ文字は綴り（`lambda`/`omega`）で代替する。R ソース自体は ASCII のはずだが `\uXXXX` エスケープで紛れ込みやすいので注意。

## ビルド・テスト手順
R は 4.6.1（`bin\x64` ではなく `bin\` 直下の実行ファイルを使う）。

```powershell
# R ファイルのみ変更後
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL .

# C++ 変更後（DLL ロック解除が必要）
Get-Process | Where-Object {$_.Name -match "Rterm"} | Stop-Process -Force
Start-Sleep -Seconds 1
Remove-Item src\*.o -Force -ErrorAction SilentlyContinue
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL .

# テスト実行（DLL コピーが必要）
Copy-Item "$env:USERPROFILE\AppData\Local\R\win-library\4.6\coresynth\libs\x64\coresynth.dll" `
  -Destination "src\coresynth.dll" -Force
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "devtools::test()"
```

修正を検証する際は、この手順でビルド・インストールしてから `devtools::test()` を実行し、DLL が更新されていることを確認すること（インストールの成否を確認しないと古い DLL でテストしてしまう）。
