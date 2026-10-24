---
trigger: glob
globs: R/*.R,src/*.cpp,DESCRIPTION,NEWS.md
---

# 手法別ドキュメントの参照先

編集・レビュー対象に対応する下記ファイルを、作業前に必ず読むこと（アルゴリズム決定の根拠・論文照合がまとまっている）。

| 対象 | 読むファイル |
|---|---|
| `R/scm.R` `src/scm.cpp` `src/optim.cpp` | `.claude/rules/scm.md` |
| `R/sdid.R` `src/sdid.cpp` | `.claude/rules/sdid.md` |
| `R/gsc.R` `src/gsc.cpp` | `.claude/rules/gsc.md` |
| `R/mc.R` `src/mc.cpp` | `.claude/rules/mc.md` |
| `R/tasc.R` `src/tasc.cpp` | `.claude/rules/tasc.md` |
| `R/si.R` `src/si.cpp` | `.claude/rules/si.md` |
| `R/scm_design.R` | `.claude/rules/scm-design.md` |
| `R/conformal.R` | `.claude/rules/conformal.md` |
| `R/scm_fit.R` `R/accessors.R` `R/broom.R` | `.claude/rules/object-model.md`（`scm_fit.R` は `scm.md` `input-validation.md` も） |
| `R/utils.R` | `.claude/rules/input-validation.md` |
| `R/plot.R` | `.claude/rules/plot.md` |
| `DESCRIPTION` `NEWS.md` | `.claude/rules/versioning.md` |
