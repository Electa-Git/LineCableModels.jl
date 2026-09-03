# Third-Party Notices

OneLineEditor 包含、使用或基于以下第三方资源。本文件汇总它们的来源、许可与使用约束。

---

## QElectroTech Elements Collection

- **来源**：QElectroTech 项目, https://qelectrotech.org
- **使用位置**：`third_party/qelectrotech/10_electric/`
- **使用方式**：作为电气元件符号的素材来源。后续会通过脚本将 `.elmt` 文件转换为 SVG + 端口元数据，集成到本项目的元件库中。
- **License**：Creative Commons Attribution 3.0 (CC-BY 3.0)
  - 全文：https://creativecommons.org/licenses/by/3.0/
  - 上游声明：见 `third_party/qelectrotech/ELEMENTS.LICENSE`

### 署名（Attribution）

> Electrical symbols originally from the QElectroTech project (https://qelectrotech.org), distributed under the Creative Commons Attribution 3.0 License.
>
> 本产品中的电气元件符号源自 QElectroTech 项目（https://qelectrotech.org），依 CC-BY 3.0 许可使用。如有改动（例如适配 Web 渲染、调整端子坐标），均为 OneLineEditor 项目所做的修改。

### 使用约束

1. **必须保留署名**：在产品的"关于"信息、文档或源代码注释中保留上述 attribution。
2. **必须声明修改**：经过转换/修改的元件 SVG 应在元数据中保留 "modified from QElectroTech original" 字样。
3. **不可暗示作者背书**：不得以任何方式暗示 QElectroTech 项目或其作者认可、推荐 OneLineEditor。
4. **禁止用作机器学习训练数据**：上游 ELEMENTS.LICENSE 显式禁止将 elements 文件用于训练 ML 模型。本项目不得、亦不得允许第三方将 `third_party/qelectrotech/` 下的 `.elmt` 文件或其衍生产物用于此目的。

### 注意

- **QElectroTech 应用程序本体**采用 GPLv2 license。本项目**不**包含或链接其应用代码，仅使用其 elements collection（CC-BY 3.0）。GPLv2 不波及本项目。
