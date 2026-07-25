---
name: build-install
description: coresynthパッケージをビルド・インストールしてテストを実行する。R/*.R またはsrc/*.cppを編集した後、変更を反映してテストを回したいときに使う。Windows特有のDLLロック・サイレントロールバックの罠を回避する正しい手順を実行する。
---

coresynth（RcppArmadillo製のRパッケージ）のビルド・インストール・テストを、
このプロジェクト固有の落とし穴を踏まずに実行する。

## 手順

変更が `R/*.R` のみの場合:

```powershell
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL .
```

変更に `src/*.cpp` が含まれる場合（DLLロック解除が必要）:

```powershell
Get-Process | Where-Object {$_.Name -match "Rterm"} | Stop-Process -Force
Start-Sleep -Seconds 1
Remove-Item src\*.o -Force -ErrorAction SilentlyContinue
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL .
```

インストール成功後、テストを実行する場合（DLLコピーが必要）:

```powershell
Copy-Item "$env:USERPROFILE\AppData\Local\R\win-library\4.6\coresynth\libs\x64\coresynth.dll" `
  -Destination "src\coresynth.dll" -Force
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "devtools::test()"
```

## 既知の落とし穴（必ず確認する）

- **`R CMD INSTALL` はサイレントロールバックする**: 失敗しても旧バージョンのDLLが
  残ったままになり、気づかずに古いビルドでテストしてしまう罠がある。
  `R CMD INSTALL .` の出力に `DONE (coresynth)` が出ているか必ず確認し、
  出ていなければテストを実行する前に原因を特定して修正する。
- **`"cannot remove earlier installation"` エラー**: 上記の `Stop-Process` 手順
  （Rtermプロセスの強制終了→1秒待機→`.o`ファイル削除）で解消する。それでも
  解消しない場合はユーザーに報告する（他のRセッションがDLLをロックしている
  可能性がある）。
- R は 4.6.1（`bin\x64` サブディレクトリではなく `bin\` 直下の実行ファイルを使用）。
