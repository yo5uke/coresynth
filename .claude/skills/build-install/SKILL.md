---
name: build-install
description: coresynthパッケージをビルド・インストールしてテストを実行する。R/*.R または src/*.cpp を編集した後、変更を反映してテストを回したいときに使う。Windows特有のDLLロック・サイレントロールバックの罠を回避する正しい手順を実行する。
---

coresynth（RcppArmadillo 製の R パッケージ）のビルド・インストール・テストを、Windows 環境固有の落とし穴を踏まずに実行する。

## 手順

R は `C:\Program Files\R\R-4.6.2\bin\` 直下の実行ファイル（`bin\x64` ではない）を使用する。

変更が `R/*.R` のみの場合:

```powershell
& "C:\Program Files\R\R-4.6.2\bin\R.exe" CMD INSTALL .
```

変更に `src/*.cpp` が含まれる場合（DLL ロック解除が必要）:

```powershell
Get-Process | Where-Object {$_.Name -match "Rterm"} | Stop-Process -Force
Start-Sleep -Seconds 1
Remove-Item src\*.o, src\*.dll -Force -ErrorAction SilentlyContinue
& "C:\Program Files\R\R-4.6.2\bin\R.exe" CMD INSTALL .
```

インストール成功後、テストを実行する場合（DLL コピーが必要）:

```powershell
Copy-Item "$env:USERPROFILE\AppData\Local\R\win-library\4.6\coresynth\libs\x64\coresynth.dll" `
  -Destination "src\coresynth.dll" -Force
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" -e "devtools::test()"
```

## 既知の落とし穴（必ず確認する）

- **`R CMD INSTALL` はサイレントロールバックする**: 失敗しても旧バージョンの DLL が残ったままになり、古いビルドでテストしてしまう罠がある。出力に `DONE (coresynth)` が出ているか必ず確認し、出ていなければ原因を修正してからテストする。
- **`DONE (coresynth)` が出ても DLL が差し替わっていないことがある**: 途中に `problem copying ...coresynth.dll ...: Permission denied` 警告が出ることがあるため、警告行も確認し、疑わしければインストール先 DLL の更新時刻を確認する。
- **DLL ロック元が `Rterm` とは限らない**: Positron の R コンソール（プロセス名 `ark`）が `coresynth` を読み込んでいるとロックされる。特定コマンド: `Get-Process | ForEach-Object { if ($_.Modules | Where-Object { $_.FileName -like "*coresynth*" }) { $_.Name } }`。ユーザー作業中のセッションは勝手に終了せず閉じてもらうよう依頼する（待てない計測は `--library=<scratch>` を使用）。
- **`devtools::test()` / `load_all()` は `src/*.o` を残す**: 次の `R CMD INSTALL .` が再利用するため、C++ 変更の有無にかかわらず速度計測前などは INSTALL 前に `Remove-Item src\*.o, src\*.dll` で統一条件を保つ。
- **`"cannot remove earlier installation"` エラー**: 上記 `Stop-Process` 手順で解消する。解消しない場合は他の R セッションがロックしている可能性があるためユーザーに報告する。
