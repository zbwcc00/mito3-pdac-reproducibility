# PDAC mitoxyperilysis evidence ledger

日期：2026-08-26

## 结论先行

当前主线不需要更换，但需要精确命名为：**PDAC 中 mitoxyperilysis-associated tumor–myeloid state 与治疗重塑、患者生存和免疫排斥的关联**。现有数据支持“状态关联链”，尚不足以支持“mitoxyperilysis 导致 ICI 耐药”的因果链，也不能把 RNA 代理分数称为真实线粒体孔化、蛋白活化或分泌活性。

## 已完成证据

| 层级 | 队列/方法 | 可用单位 | 主要观察 | 证据定位 |
|---|---|---:|---|---|
| 免疫治疗 RNA | PRINCE/Padrón 2022 | 65 RNA，53 个应答，63 个 OS | nivolumab-only 中线粒体氧化模块 Δmedian≈−0.499，Wilcoxon P≈0.0186，FDR≈0.065；全队列 OS HR≈1.61，P≈0.040，FDR≈0.278 | 重要名义支持，未达稳定 ICI 验证 |
| 免疫治疗外周蛋白 | PRINCE Olink | 103 基线血浆 | 应答比较所有蛋白 FDR≈0.981；OS 层 IL-6、IL-8、HGF、MMP12、MMP7、CCL7、VEGF-A 等有 FDR 信号 | 外周血辅助层，不能替代肿瘤组织蛋白/空间 |
| 新辅助治疗 Nanostring | GSE129492 | 24 例，4 组各 6 例 | 线粒体氧化、先天输入、炎症背景和髓系生态模块在治疗组间改变；XRT 对线粒体氧化、先天输入及炎症背景的差异经 BH 校正仍存在 | 治疗重塑支持，不是 ICI 反应验证 |
| 空间转录组 | GSE240078 | 223 AOI | GEO/论文明确 36 例、4 对纵向样本；GEO 标题前缀解析出 40 个候选 ID，尚未获得可靠 36 例映射 | 空间方向可用；患者级结果暂不纳入最终独立验证 |
| 肿瘤细胞遗传依赖 | DepMap 24Q4 | 64 个 PDAC 模型被标注，41 个目标列完整 | RICTOR 中位依赖概率≈0.459，17/41≥0.5；RHOA≈0.423，17/41≥0.5；BAK1≈0.370，13/41≥0.5；BAX/BID/MYD88/SPP1/IL1B 均低 | 候选依赖性排序；不是 TME 或 ICI 证据 |
| 肿瘤细胞药敏 | PRISM 19Q4 | 37 个 PDAC 细胞系，4686 个处理列 | mTOR/PI3K 类药物敏感性高度异质；BGT226、GSK2126458 等位于全药物敏感性前列；BAX 激活剂 BAM7 居中；ROCK 抑制剂分散 | 药理学辅助，不能归因于 RICTOR 或 mitoxyperilysis |
| 单细胞虚拟扰动 | GSE154778 | 1160 个 singlet 髓系细胞，15 个患者 | 线性模型校正文库量、髓系亚群和患者后，RHOA/BAX 的预测模块改变较明显；RICTOR 检出率仅约 13%，多数模拟 KO 的中位变化为 0 | 透明敏感性筛选；不是 CellOracle/scTenifoldKnk 因果结果 |

## 主线审查

### 保留部分

1. **PDAC 肿瘤—髓系生态—生存**仍是最稳健的主线。髓系定位、空间/区域分层和患者级生存分析可以组成正文主体。
2. **mitoxyperilysis**适合作为机制启发和状态标签，而不是已经被 PDAC 临床样本直接证明的死亡方式。
3. **ICI**应放在转化假设和外部支持层。PRINCE 的 RNA 结果有方向，但 FDR 和跨队列重复尚未达到“稳定预测标志物”标准。
4. **DepMap/PRISM/虚拟扰动**适合做候选优先级和机制三角验证，不能替代蛋白、空间共定位或干预实验。

### 必须降级或避免的表述

- 不写“mitoxyperilysis causes ICI resistance”。
- 不写“RICTOR is a validated PDAC dependency”；当前只有约 41% 的可用模型达到经验阈值，且未做表达/谱系分层与外部复现。
- 不把 Olink 外周血蛋白写成肿瘤局部 SecAct 分泌活性。
- 不把当前透明线性虚拟扰动命名为 CellOracle/scTenifoldKnk 结果；若加入正文，应称为 model-based in-silico perturbation sensitivity analysis。
- 不把 GSE240078 的 40 个标题前缀当作 40 个患者。论文设计是 36 例，必须得到补充表或作者映射后才能做患者级独立验证。

## 是否继续加入 CellOracle/scTenifoldKnk

目前不建议把完整 CellOracle/scTenifoldKnk 作为主结果。原因是 RICTOR 在 GSE154778 髓系细胞中检出率低，现有数据没有匹配的 scATAC/TF motif 先验，且虚拟扰动会高度依赖网络先验。当前透明模型已足够用于候选排序；只有在获得可靠细胞类型特异 GRN、足够目标基因检出率或外部 Perturb-seq 证据后，才值得升级为正式工具分析。

## LINCS 是否必需

暂不纳入主分析。PRISM 已提供 PDAC 细胞系药敏层；LINCS 的广泛转录响应签名缺少本项目患者/细胞背景匹配，加入后更可能增加机制解释空间而非增加独立验证强度。若后续需要药物重定位图，可将 LINCS 作为补充图，并明确其为跨细胞背景的转录响应推断。

## 六步执行状态

1. GSE129492 Nanostring：完成。
2. GSE240078 患者映射边界与 PRINCE Olink FDR：完成，空间患者映射保留为阻断项。
3. DepMap 24Q4 依赖性：完成。
4. PRISM 19Q4 药敏：完成；LINCS 暂不做，原因见上。
5. 单细胞虚拟扰动可行性：完成，采用透明模型而非冒充 CellOracle/scTenifoldKnk。
6. 证据账本与主线更新：完成。

## 下一阶段优先级

1. 取得 GSE240078 36 例临床/ROI 映射后，重新做患者级空间验证。
2. 对 PRINCE RNA、Olink 和现有单细胞结果做预注册式方向一致性汇总，而不是继续堆叠名义 P 值。
3. 若能获得 CPTAC 配对磷酸化或组织蛋白，优先验证 pAKT-S473/mTORC2、BAX/BAK1 和髓系空间共定位。
4. 将 DepMap/PRISM/虚拟扰动放入机制补充层，候选优先级暂定为 RICTOR、RHOA、BAK1，而不是宣称已验证驱动因子。
