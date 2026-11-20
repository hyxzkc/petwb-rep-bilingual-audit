# -*- coding: utf-8 -*-
"""
organ_state_extractor.py  (v2.0 — 正式分析版)
=============================================
基于 terminology_ontology_v01.md 的中英 PET/CT 报告「19 器官 x 4 状态 + 摄取命题」规则抽取器。

方法学口径（与 GPT 讨论后锁定，逐条落实于代码）：
  修正1  四状态主标签之外，新增 uptake_proposition（摄取命题）层。
         中文"摄取未见增高/未见异常增高"= not_increased（主状态仍 Normal），
         "摄取未见异常"= not_abnormal，"生理性摄取"= physiological，
         明确"摄取正常"= normal_explicit；英文同理。
         新增语义漂移 flag = over_normalization：
           zh proposition ∈ {not_increased, not_abnormal, physiological}
           且 en proposition == normal_explicit
           （中文只否定异常，英文却断言"完全正常"，语义特异性损失）。
  修正2  区分四种 Unknown 成因：not_mentioned / unknown_due_to_extraction_failure /
         uncertain_in_report；已知状态则该列留空。四状态主标签不变。
  修正3  分层抽取结构，五步显式化：
         ① 器官定位 ② 事件识别 ③ 否定/不确定范围 ④ 程度与数值 ⑤ 四状态映射 + 语义限定映射。
         每个状态单元保留 evidence（中英片段）、trigger（触发词）、
         rule_version（规则版本）、confidence（置信度）。
  修正4  器官覆盖三层审计（organ_coverage_audit.json）：
         ① 原始词法覆盖差异 ② 同义词归一化后覆盖差异 ③ 正则标记需人工确认的差异。
         只有 ③ 可能进入漂移结论。

设计原则（方法学预实验→正式分析）：
1. 只读真实报告，不补写、不推测；抽不出一律 Unknown，evidence 留空或保留段落供人工。
2. 四状态 = Increased / Decreased / Absent / Normal；Unknown 独立保留，不自动归一为 Normal。
3. 抽取成功率与两个一致率口径固定，绝不混用：
   - 可比口径：两侧均非 Unknown 的器官对；
   - 含 both_unknown 口径：两侧均 Unknown 也计为一致。
4. 文本规范化：_x000D_ / NBSP / 全角连字符。

用法：
    python organ_state_extractor.py            # 30 对验证 → state_alignment30.csv / p0_pilot_stats.json
    python organ_state_extractor.py --full     # 490 对 → state_alignment490.csv / full_stats.json /
                                               #          organ_coverage_audit.json
"""
import argparse
import csv
import json
import os
import re
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "petwb_data")
RULE_VERSION = "v2.0_20261120"

# ---------------------------------------------------------------------------
# 一、19 器官：中英锚点（regex）。顺序与本体一致。
#     锚点只用于"定位包含该器官的句子"（修正3 第①步），不参与状态判定。
# ---------------------------------------------------------------------------
ORGAN_SPECS = [
    ("01 Cranium and Brain (颅脑)",
     [r"脑", r"颅脑", r"大脑", r"小脑", r"脑干", r"脑室", r"脑沟", r"脑裂", r"脑池",
      r"脑实质", r"脑白质", r"基底节", r"鞍区", r"四叠体", r"胼胝体", r"松果体"],
     [r"\bbrain\b", r"\bcerebr\w*\b", r"\bcranial\b", r"\bventricles?\b",
      r"\bcerebell\w*\b", r"\bpons\b", r"\bbasal ganglia\b", r"\bwhite matter\b"]),

    ("02 Eyeballs (眼球)",
     [r"眼球", r"眼眶", r"球后", r"视神经", r"双眼"],
     [r"\beyeballs?\b", r"\borbit\w*\b", r"\bretrobulbar\b", r"\boptic nerves?\b", r"\beyes?\b"]),

    ("03 Nasal Cavity and Sinuses (鼻腔与鼻窦)",
     [r"副鼻窦", r"鼻窦", r"上颌窦", r"筛窦", r"蝶窦", r"额窦", r"鼻腔", r"鼻中隔", r"鼻甲"],
     [r"\bsinuses?\b", r"\bmaxillary\b", r"\bethmoid\b", r"\bsphenoid\b",
      r"\bfrontal sinus\b", r"\bnasal cavity\b"]),

    ("04 Pharynx and Parapharyngeal Space (咽部及咽旁间隙)",
     [r"鼻咽", r"口咽", r"喉咽", r"咽隐窝", r"咽旁间隙", r"咽鼓管", r"咽部", r"咽壁"],
     [r"\bnasopharynx\b", r"\boropharynx\b", r"\blaryngopharynx\b", r"\bpharyngeal\b",
      r"\bparapharyngeal\b", r"\bEustachian\b", r"\bpharynx\b"]),

    ("05 Palatine Tonsils and Larynx (腭扁桃体与喉)",
     [r"扁桃体", r"腭扁桃体", r"喉部", r"声门", r"声带", r"喉咽"],
     [r"\btonsils?\b", r"\blarynx\b", r"\bglottis\b", r"\bvocal cords?\b"]),

    ("06 Salivary Glands and Thyroid (唾液腺与甲状腺)",
     [r"腮腺", r"颌下腺", r"舌下腺", r"唾液腺", r"甲状腺"],
     [r"\bparotid\b", r"\bsubmandibular glands?\b", r"\bsublingual\b",
      r"\bsalivary\b", r"\bthyroid\b"]),

    ("07 Cervical Lymph Nodes (颈部淋巴结)",
     [r"颈深间隙", r"颈深上间隙", r"颈部淋巴", r"颈淋巴结", r"颈后三角",
      r"颏下", r"颌下及颏下", r"颈部间隙"],
     [r"\bdeep cervical\b", r"\bcervical (spaces?|lymph|nodes?)\b",
      r"\bposterior cervical triangle\b", r"\bsubmental\b", r"\bneck lymph\b",
      r"\bsubmandibular region\b"]),

    ("08 Lungs and Thoracic Cavity (肺与胸腔)",
     [r"肺", r"胸膜", r"胸腔", r"肺门", r"肺野", r"支气管", r"肺实质", r"肺叶", r"肺结节"],
     [r"\blungs?\b", r"\bpleura\w*\b", r"\bthoracic cavity\b", r"\bpulmonary\b",
      r"\bbronch\w*\b", r"\bhila\b", r"\bhilum\b", r"\bpneumothorax\b"]),

    ("09 Mediastinum and Heart (纵隔与心脏)",
     [r"纵隔", r"心影", r"心脏", r"心肌", r"心包", r"主动脉", r"大血管", r"纵隔血池"],
     [r"\bmediastin\w*\b", r"\bheart\b", r"\bcardiac\b", r"\bpericard\w*\b",
      r"\baorta\w*\b", r"\bcardiac silhouette\b"]),

    ("10 Esophagus (食管)",
     [r"食管(?!旁|右旁|左旁|气管)", r"食道"],
     [r"\besophagus\b", r"\boesophagus\b", r"\besophageal\b", r"\boesophageal\b"]),

    ("11 Liver (肝脏)",
     [r"肝脏", r"肝实质", r"肝内", r"肝缘", r"肝裂", r"肝右叶", r"肝左叶",
      r"肝门部", r"肝门区", r"肝实质内"],
     [r"\bliver\b", r"\bhepatic\b"]),

    ("12 Gallbladder (胆囊)",
     [r"胆囊", r"胆道", r"胆总管"],
     [r"\bgallbladder\b", r"\bbiliary tract\b", r"\bcommon bile duct\b"]),

    ("13 Pancreas (胰腺)",
     [r"胰腺", r"胰头", r"胰体", r"胰尾", r"胰管", r"主胰管"],
     [r"\bpancreas\b", r"\bpancreatic\b"]),

    ("14 Spleen (脾脏)",
     [r"脾脏", r"脾旁", r"脾实质", r"脾(?!脏)"],
     [r"\bspleen\b", r"\bsplenic\b"]),

    ("15 Kidneys and Adrenal Glands (肾脏与肾上腺)",
     [r"肾脏", r"双肾", r"左肾", r"右肾", r"肾上腺", r"肾盂", r"肾盏", r"输尿管", r"肾实质"],
     [r"\bkidneys?\b", r"\brenal\b", r"\badrenals?\b", r"\bureters?\b"]),

    ("16 Gastrointestinal Tract (胃肠道)",
     [r"(?<!肝)胃", r"肠道", r"小肠", r"结肠", r"直肠", r"肠管", r"消化道", r"十二指肠", r"肛管"],
     [r"\bstomach\b", r"\bgastric\b", r"\bbowel\b", r"\bintestin\w*\b", r"\bcolon\b",
      r"\brectum\b", r"\bgastrointestinal\b", r"\bduoden\w*\b", r"\banal canal\b"]),

    ("17 Prostate/Uterus and Bladder (前列腺/子宫与膀胱)",
     [r"前列腺", r"子宫", r"膀胱", r"附件区", r"双侧附件", r"单侧附件", r"卵巢"],
     [r"\bprostate\b", r"\buterus\b", r"\bbladder\b", r"\badnexa\b", r"\bovaries\b"]),

    ("18 Abdominal and Pelvic Cavities (腹腔与盆腔)",
     [r"腹腔", r"盆腔", r"腹膜后", r"腹膜", r"腹水", r"腹股沟"],
     [r"\babdominal cavity\b", r"\babdomen\b", r"(?<!renal )\bpelvis\b", r"\bpelvic\b",
      r"\bperitoneum\b", r"\bperitoneal\b", r"\bretroperitone\w*\b"]),

    ("19 Spine and Bones (脊柱与骨骼)",
     [r"脊柱", r"脊椎", r"椎体", r"椎间盘", r"骨骼", r"骨盆", r"髂骨", r"肋骨", r"胸骨",
      r"肩胛骨", r"骨转移", r"全身骨骼"],
     [r"\bspine\b", r"\bspinal\b", r"\bvertebra\w*\b", r"\bbones?\b", r"\bskeleton\b",
      r"\bribs?\b", r"\bsternum\b", r"\bscapula\b", r"\bsacrum\b"]),
]

# 修正4 用：术语本体 v0.1 第一节的同义词（归一化后覆盖比较）
ONTOLOGY_SYNONYMS = {
    1:  (["颅骨与脑", "头颅", "脑部", "颅内"], ["skull and brain", "head", "cranium"]),
    2:  (["双眼", "眼眶内眼球"], ["eyes", "orbits"]),
    3:  (["鼻腔及副鼻窦", "鼻旁窦", "副鼻窦"], ["nasal cavity", "paranasal sinuses"]),
    4:  (["鼻咽", "口咽", "喉咽", "咽旁间隙"], ["nasopharynx", "oropharynx", "hypopharynx"]),
    5:  (["扁桃体", "喉部", "声门"], ["tonsils", "larynx", "glottis"]),
    6:  (["腮腺", "颌下腺", "舌下腺", "甲状腺"], ["parotid", "submandibular gland", "thyroid"]),
    7:  (["颈淋巴结", "颈部淋巴引流区"], ["neck lymph nodes", "cervical LN"]),
    8:  (["双肺", "肺野", "胸膜腔", "胸腔"], ["lungs", "thoracic cavity", "pleura"]),
    9:  (["纵隔", "心影", "心肌", "心包"], ["mediastinum", "heart", "pericardium"]),
    10: (["食道"], ["oesophagus"]),
    11: (["肝", "肝实质"], ["hepatic"]),
    12: (["胆囊窝", "胆道"], ["gallbladder", "biliary tract"]),
    13: (["胰", "胰头", "胰体", "胰尾"], ["pancreatic"]),
    14: (["脾"], ["splenic"]),
    15: (["双肾", "肾", "肾上腺"], ["kidneys", "adrenals"]),
    16: (["胃", "小肠", "结肠", "直肠", "消化道"], ["GI tract", "stomach", "bowel", "colon"]),
    17: (["前列腺", "子宫", "附件", "膀胱"], ["prostate", "uterus", "bladder", "adnexa"]),
    18: (["腹腔", "盆腔", "腹膜", "腹膜后"], ["abdominal cavity", "pelvis", "peritoneum", "retroperitoneum"]),
    19: (["脊椎", "椎体", "骨盆骨", "骨转移"], ["spine", "vertebrae", "skeleton", "bone"]),
}

# 修正4 第③层：这些词使器官锚点归属不确定（跨器官泄漏/解剖边界歧义），需人工确认
AMBIGUOUS_MARKERS = {
    "zh": [r"旁", r"间隙", r"窝", r"淋巴结", r"骨髓", r"肩周", r"颅底", r"膈", r"附件", r"切缘"],
    "en": [r"lymph node", r"renal pelvis", r"adnexal", r"bone marrow", r"skull base",
           r"axillary", r"hilar", r"parotid", r"submandibular"],
}

# ---------------------------------------------------------------------------
# 二、摄取命题（uptake proposition）词表 —— 修正1 核心
#     每条 = (regex, proposition)。列表顺序即命中优先级：否定/生理/正常在前，
#     异常（increased/decreased/absent）在后，保证"未见增高"不会被判为 Increased。
#     proposition → 四状态主标签的映射见 PROP_TO_STATE。
# ---------------------------------------------------------------------------
PROP_TO_STATE = {
    "increased":       "Increased",
    "decreased":       "Decreased",
    "absent":          "Absent",
    "normal_explicit": "Normal",
    "not_increased":   "Normal",
    "not_abnormal":    "Normal",
    "physiological":   "Normal",
    "unknown":         "Unknown",
}
PROPOSITIONS = ["increased", "decreased", "absent", "normal_explicit",
                "not_increased", "not_abnormal", "physiological", "unknown"]

ZH_PROP_RULES = [
    # --- 否定"异常增高"：not_increased（必须先于"未见异常"） ---
    (r"摄取(均|亦)?未见(明显|局灶性)?异常(增高|升高|增加|浓聚)", "not_increased"),
    (r"代谢(均|亦)?未见(明显|局灶性)?异常(增高|升高|增加)", "not_increased"),
    (r"摄取(均|亦)?未见(明显|局灶性)?(增高|升高|增加|浓聚)", "not_increased"),
    (r"代谢(均|亦)?未见(明显|局灶性)?(增高|升高|增加)", "not_increased"),
    (r"未见(明显|局灶性)?异常(增高|升高|增加)的?(FDG)?(摄取|代谢)", "not_increased"),
    (r"(摄取|代谢)(不高|未见增高)", "not_increased"),
    # --- 否定异常：not_abnormal ---
    (r"(摄取|代谢)(均|亦)?未见(明显|局灶性)?异常", "not_abnormal"),
    (r"(未见|未发现)(明显)?异常(FDG)?(摄取|代谢)", "not_abnormal"),
    (r"(无|没有)(明显)?异常(FDG)?(摄取|代谢)", "not_abnormal"),
    (r"FDG(摄取|代谢)?无异常", "not_abnormal"),
    (r"无异常(FDG)?(摄取|代谢)", "not_abnormal"),
    (r"(摄取|代谢)未见异常", "not_abnormal"),
    # --- 生理性 / 本底：physiological ---
    (r"生理性(FDG)?(摄取|代谢|浓聚)", "physiological"),
    (r"(FDG)?(摄取|代谢)(呈|为)? ?(生理性|本底|基础)摄取?", "physiological"),
    (r"本底摄取", "physiological"),
    # --- 明确正常：normal_explicit ---
    (r"摄取正常", "normal_explicit"),
    (r"代谢正常", "normal_explicit"),
    # --- 异常：increased / decreased / absent ---
    (r"摄取缺失", "absent"), (r"摄取消失", "absent"), (r"无(FDG)?摄取", "absent"),
    (r"未见(FDG)?摄取", "absent"), (r"无(示踪剂|放射性)摄取", "absent"),
    (r"FDG摄取缺如", "absent"), (r"(FDG)?代谢缺如", "absent"),
    (r"摄取(较前)?(减低|降低|减少|下降)", "decreased"),
    (r"代谢(较前)?(减低|降低|减少|下降)", "decreased"),
    (r"低摄取", "decreased"),
    (r"摄取(轻度|稍|轻微|欠均匀|不均匀|明显|显著|局灶性)?(增高|升高|增加|活跃)", "increased"),
    (r"代谢(轻度|稍|轻微|欠均匀|不均匀|明显|显著|局灶性)?(增高|升高|增加)", "increased"),
    (r"(轻度|稍|轻微)FDG(摄取|代谢)", "increased"),
    (r"高摄取", "increased"), (r"高代谢", "increased"), (r"轻微摄取", "increased"),
    (r"异常浓聚", "increased"), (r"浓聚", "increased"), (r"FDG高摄取", "increased"),
]

EN_PROP_RULES = [
    # --- not_increased ---
    (r"\bno\s+(abnormally|significantly|markedly)?\s*(increased|elevated|high)\s+(FDG\s+)?(uptake|metabolism)\b", "not_increased"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(was|is|were|are|showed|shows)?\s*not\s+(abnormally|significantly|markedly)?\s*(increased|elevated|high)\b", "not_increased"),
    (r"\bno\s+(significant\s+|focal\s+)?(abnormal\s+)?increase\s+in\s+(FDG\s+)?(uptake|metabolism)\b", "not_increased"),
    (r"\b(FDG\s+)?(uptake|metabolism)\b.{0,50}?\bnot\s+abnormally\s+increased\b", "not_increased"),
    # --- not_abnormal ---
    (r"\bno\s+(obvious|significant|obviously|significantly)?\s*abnormal\s+(FDG\s+)?(uptake|metabolism)\b", "not_abnormal"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(was|is|were|are|showed|shows)?\s*(not|no longer)\s+(significantly\s+)?abnormal\b", "not_abnormal"),
    (r"\b(FDG\s+)?uptake\s+was\s+not\s+abnormal\b", "not_abnormal"),
    (r"\bno\s+(abnormal|significant|high|elevated|obvious)\s+(local\s+)?(FDG\s+)?(uptake|metabolism)\b", "not_abnormal"),
    (r"\bwithout\s+(abnormal|increased|significant)\s+(FDG\s+)?(uptake|metabolism)\b", "not_abnormal"),
    (r"\bno\s+(significant\s+|obvious\s+)?abnormalit(y|ies)\b.{0,80}?\b(FDG\s+)?(uptake|metabolism)\b", "not_abnormal"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(showed|shows)\s+no\s+(significant\s+)?abnormalit(y|ies|increase)\b", "not_abnormal"),
    (r"\bno\s+abnormal\b.{0,40}?\b(FDG\s+)?(uptake|metabolism)\b", "not_abnormal"),
    # --- physiological ---
    (r"\bphysiological\s+(FDG\s+)?(uptake|metabolism|accumulation)\b", "physiological"),
    (r"\b(FDG\s+)?(uptake|metabolism)\b.{0,35}?\bphysiological\b", "physiological"),
    (r"\bbackground\s+(FDG\s+)?(uptake|metabolism)\b", "physiological"),
    # --- normal_explicit ---
    (r"\bnormal\s+(FDG\s+)?(uptake|metabolism|accumulation)\b", "normal_explicit"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(was|is|were|are|showed|shows)\s+normal\b", "normal_explicit"),
    (r"\b(FDG\s+)?(uptake|metabolism)\b.{0,40}?\b(is|was|were|are)\s*normal\b", "normal_explicit"),
    (r"\bnormal\s+in\b.{0,60}?\b(FDG\s+)?(uptake|metabolism)\b", "normal_explicit"),
    # --- absent ---
    (r"\bno\s+FDG\s+uptake\b", "absent"),
    (r"\babsent\s+(FDG\s+)?(uptake|metabolism)\b", "absent"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(is|was|are|were)\s+absent\b", "absent"),
    (r"\bno\s+(tracer|radioactivity)\s+uptake\b", "absent"),
    # --- decreased ---
    (r"\b(reduced|decreased|diminished|low|lower)\s+(FDG\s+)?(uptake|metabolism)\b", "decreased"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(was|is|were|are)?\s*(reduced|decreased|low)\b", "decreased"),
    # --- increased ---
    (r"\b(increased|elevated|high|heightened|intense|marked|hypermetabolic)\s+(FDG\s+)?(uptake|metabolism)\b", "increased"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(was|is|were|are)?\s*(unevenly|slightly|mildly|markedly|significantly|abnormally)?\s*(increased|elevated|high|marked)\b", "increased"),
    (r"\b(mild|slight|moderate)\s+(FDG\s+)?(uptake|metabolism)\b", "increased"),
    (r"\b(FDG\s+)?(uptake|metabolism)\s+(elevation|elevations)\b", "increased"),
    (r"\bFDG-high\s+uptake\b", "increased"), (r"\bincrease\s+in\s+FDG\s+uptake\b", "increased"),
    (r"\b(FDG|tracer)\s+accumulation\b", "increased"), (r"\bhypermetabolic\b", "increased"),
    (r"\bshowed\s+(increased|elevated|high)\s+FDG\b", "increased"),
]

def _compile_prop(pairs):
    return [(re.compile(p, re.IGNORECASE), prop) for p, prop in pairs]

PROP_RULES = {
    "zh": _compile_prop(ZH_PROP_RULES),
    "en": _compile_prop(EN_PROP_RULES),
}

# 两级优先级，用途不同，切勿混用：
# 1) 子句级：Normal 状态优先 —— 解析"摄取未见增高"（=Normal + not_increased，不是 Increased）
CLAUSE_PRIORITY = ["Normal", "Absent", "Decreased", "Increased"]
# 2) 器官级：Increased 优先 —— 器官内出现高代谢病灶时取 dominant abnormal finding
ORGAN_PRIORITY = ["Increased", "Decreased", "Absent", "Normal"]

# ---------------------------------------------------------------------------
# 三、临床限定词（修正3 第④步：程度与数值 / 第⑤步：语义限定映射，单独管理）
# ---------------------------------------------------------------------------
QUALIFIERS = {
    "negation":    _compile_prop([(r"未见", "x"), (r"未发现", "x"), (r"无", "x"), (r"阴性", "x"),
                                  (r"排除", "x"), (r"没有", "x"), (r"\bno\b", "x"),
                                  (r"\bnone\b", "x"), (r"\bwithout\b", "x"),
                                  (r"\bnegative\b", "x"), (r"rule out", "x")]),
    "uncertainty": _compile_prop([(r"可疑", "x"), (r"可能", "x"), (r"考虑", "x"), (r"待排", "x"),
                                  (r"不除外", "x"), (r"倾向于", "x"), (r"\bsuspicious\b", "x"),
                                  (r"\bpossible\b", "x"), (r"\blikely\b", "x"),
                                  (r"\bconsider\b", "x"), (r"cannot exclude", "x"),
                                  (r"\bprobably\b", "x")]),
    "degree":      _compile_prop([(r"轻度", "x"), (r"中度", "x"), (r"重度", "x"), (r"明显", "x"),
                                  (r"少许", "x"), (r"稍", "x"), (r"\bmild\b", "x"),
                                  (r"\bmoderate\b", "x"), (r"\bsevere\b", "x"),
                                  (r"\bslightly\b", "x"), (r"\bmarked\b", "x"),
                                  (r"\bsignificant\b", "x")]),
    "laterality":  _compile_prop([(r"左侧", "x"), (r"右侧", "x"), (r"双侧", "x"), (r"单发", "x"),
                                  (r"\bleft\b", "x"), (r"\bright\b", "x"),
                                  (r"\bbilateral\b", "x"), (r"\bunilateral\b", "x")]),
    "numeric":     _compile_prop([(r"SUVmax\s*=?\s*\d", "x"), (r"SUVpeak\s*=?\s*\d", "x"),
                                  (r"\d+(\.\d+)?\s*(mm|cm)", "x"),
                                  (r"\d+\s*[×*x]\s*\d+", "x")]),
}
UNCERTAINTY_RX = [rx for rx, _ in QUALIFIERS["uncertainty"]]

# ---------------------------------------------------------------------------
# 四、工具函数（修正3 第①步前置：文本规范化）
# ---------------------------------------------------------------------------
def normalize(text):
    """统一换行标记、NBSP、全角/特殊连字符、空白。"""
    if text is None:
        return ""
    text = text.replace("_x000D_", "\n").replace("\r\n", "\n").replace("\r", "\n")
    text = text.replace("\xa0", " ")
    for ch in ("\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\uff0d", "\u2212"):
        text = text.replace(ch, "-")
    text = re.sub(r"[ \t]+", " ", text)
    return text

def split_sentences(text, lang):
    if lang == "zh":
        parts = re.split(r"[。；;\n]+", text)
    else:
        parts = re.split(r"(?<=[.;])\s+|[\n]+", text)
    return [p.strip() for p in parts if p and p.strip()]

def split_subclauses(sentence, lang):
    """把长句再切成逗号级子句，避免同一句内多个器官/多种状态互相覆盖。"""
    sep = r"[，,]" if lang == "zh" else r","
    return [p.strip() for p in re.split(sep, sentence) if p and p.strip()]

def classify_clause(clause, lang):
    """修正3 第②③⑤步：返回该子句的单一最佳 (state, proposition, trigger)。

    同一否定短语（如 "no increased FDG uptake"）会同时命中 not_increased 与
    increased 两条规则，必须只取子句级优先级最高的一条（Normal 命题在前），
    否则器官级 Increased 优先会把否定式错误地升级为 Increased。
    """
    hits = []
    for rx, prop in PROP_RULES[lang]:
        m = rx.search(clause)
        if m:
            hits.append((PROP_TO_STATE[prop], prop, m.group(0).strip()))
    if not hits:
        return []
    order = {s: i for i, s in enumerate(CLAUSE_PRIORITY)}
    hits.sort(key=lambda x: order.get(x[0], 9))
    return [hits[0]]

def extract_qualifiers(clauses, lang):
    """从子句中抽取临床限定词，返回 {类别: [去重后的命中词]}。"""
    blob = " ".join(clauses)
    out = {}
    for cat, pairs in QUALIFIERS.items():
        hits = []
        for rx, _ in pairs:
            m = rx.search(blob)
            if m:
                tok = m.group(0).strip()
                if tok not in hits:
                    hits.append(tok)
        out[cat] = hits
    return out

def qualifiers_to_str(q):
    return "; ".join("{}:{}".format(k, "/".join(v)) for k, v in q.items() if v)

def snip(text, limit=160):
    t = re.sub(r"\s+", " ", text).strip()
    return t if len(t) <= limit else t[:limit] + "…"

def decide_unknown_reason(has_clauses, clauses, lang):
    """修正2：四状态主标签为 Unknown 时，区分三种成因。

    - 段落关键词不存在           -> not_mentioned
    - 段落存在但状态词未命中     -> unknown_due_to_extraction_failure（保留 evidence）
    - 有描述但用不确定表达       -> uncertain_in_report
    """
    if not has_clauses:
        return "not_mentioned"
    blob = " ".join(clauses)
    if any(rx.search(blob) for rx in UNCERTAINTY_RX):
        return "uncertain_in_report"
    return "unknown_due_to_extraction_failure"

def decide_confidence(state, states_seen, unknown_reason, qualifiers):
    """置信度（启发式，0~1，保留两位）：清晰单元高、混合/不确定单元低。"""
    if state == "Unknown":
        return {"not_mentioned": "0.90",
                "unknown_due_to_extraction_failure": "0.20",
                "uncertain_in_report": "0.40"}.get(unknown_reason, "0.20")
    conf = 0.95
    if len(states_seen) > 1:
        conf -= 0.20                       # 器官内多种状态并存，取 dominant
    if qualifiers.get("uncertainty"):
        conf -= 0.15                       # 有 hedge 描述
    return "{:.2f}".format(max(conf, 0.30))

# ---------------------------------------------------------------------------
# 五、核心：单篇报告 → {器官: 状态单元}（修正3：五步分层抽取）
# ---------------------------------------------------------------------------
def extract_report(text, lang):
    # ① 器官定位：锚点 regex + 无锚点续句窗口（英文报告常把器官与摄取短语分号拆开）
    text = normalize(text)
    sentences = split_sentences(text, lang)
    anchors = [[re.compile(p, re.IGNORECASE) for p in (zhs if lang == "zh" else ens)]
               for _, zhs, ens in ORGAN_SPECS]
    sent_organs = []
    for s in sentences:
        hit = set()
        for i, rxs in enumerate(anchors):
            if any(rx.search(s) for rx in rxs):
                hit.add(i)
        sent_organs.append(hit)

    result = {}
    for idx, (name, _, _) in enumerate(ORGAN_SPECS):
        clauses, seen = [], set()
        for j, s in enumerate(sentences):
            if idx not in sent_organs[j]:
                continue
            candidates = [s]
            if j + 1 < len(sentences) and not sent_organs[j + 1]:
                candidates.append(sentences[j + 1])
            for c in candidates:
                key = c[:50]
                if key not in seen:
                    seen.add(key)
                    clauses.append(c)
        has_clauses = bool(clauses)

        # ② 事件识别 + ③ 否定/不确定范围 + ⑤ 状态映射：逗号级子句逐条判定
        #    每状态保留该器官内全部 (state, proposition, trigger, clause) 供 dominant 选择
        per_state = defaultdict(list)
        states_seen = set()
        for c in clauses:
            for sub in split_subclauses(c, lang):
                for state, prop, trig in classify_clause(sub, lang):
                    states_seen.add(state)
                    per_state[state].append({"prop": prop, "trigger": trig, "clause": c})

        state = None
        for cand in ORGAN_PRIORITY:          # 器官级：dominant abnormal finding 优先
            if per_state.get(cand):
                state = cand
                break
        is_known = state is not None

        # ④ 程度与数值 + ⑤ 语义限定映射（独立列，不混入状态）
        qualifiers = extract_qualifiers(clauses, lang)

        if is_known:
            chosen = per_state[state][0]     # 同状态多命中：取子句级优先级最高者
            proposition = chosen["prop"]
            trigger = chosen["trigger"]
            evidence = []
            for unit in per_state[state][:2]:
                evidence.append(snip(unit["clause"]))
            unknown_reason = ""
        else:
            proposition, trigger, evidence = "unknown", "", []
            unknown_reason = decide_unknown_reason(has_clauses, clauses, lang)
            # 修正2：段落存在但状态词未命中时保留 evidence 供人工
            if unknown_reason in ("unknown_due_to_extraction_failure", "uncertain_in_report"):
                evidence = [snip(c) for c in clauses[:2]]

        result[name] = {
            "state": state if is_known else "Unknown",
            "proposition": proposition,
            "trigger": trigger,
            "evidence": " | ".join(evidence),
            "qualifiers": qualifiers,
            "states_seen": sorted(states_seen),
            "unknown_reason": unknown_reason,
            "confidence": decide_confidence(
                state if is_known else "Unknown", states_seen, unknown_reason, qualifiers),
            "rule_version": RULE_VERSION,
        }
    return result

# ---------------------------------------------------------------------------
# 六、对齐与 flag 判定（修正1：新增 over_normalization）
# ---------------------------------------------------------------------------
KNOWN = {"Increased", "Decreased", "Absent", "Normal"}
OVER_NORMALIZATION_ZH = {"not_increased", "not_abnormal", "physiological"}

def is_over_normalization(z, e):
    return (z["proposition"] in OVER_NORMALIZATION_ZH
            and e["proposition"] == "normal_explicit")

def decide_state_flag(zh, en):
    """四状态主标签对齐 flag（不含语义命题层）。"""
    zs, es = zh["state"], en["state"]
    if zs == "Unknown" and es == "Unknown":
        return "both_unknown"
    if zs == "Unknown":
        return "zh_unknown"
    if es == "Unknown":
        return "en_unknown"
    if zs == es:
        return "match"
    zq = zh["qualifiers"]; eq = en["qualifiers"]
    if zq.get("uncertainty") or eq.get("uncertainty"):
        return "uncertainty_shift"
    if ({zs, es} & {"Absent", "Normal"}) and ({zs, es} & {"Increased", "Decreased"}):
        return "negation_mismatch"
    if {zs, es} == {"Increased", "Decreased"}:
        return "state_reversal"
    if zq.get("degree") != eq.get("degree") and (zq.get("degree") or eq.get("degree")):
        return "severity_shift"
    return "state_mismatch"

# ---------------------------------------------------------------------------
# 七、主流程
# ---------------------------------------------------------------------------
def load_reports():
    zh, en = {}, {}
    with open(os.path.join(DATA, "report_zh_v3.csv"), encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            zh[row["影像号"]] = row["影像描述"]
    with open(os.path.join(DATA, "report_en.csv"), encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            en[row["Image ID"]] = row["Image description"]
    return zh, en

def load_meta():
    meta = {}
    with open(os.path.join(DATA, "meta_data.csv"), encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            meta[row["Image ID"]] = row["Disease"]
    return meta

def load_sample(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def run(sample, zh_reports, en_reports):
    rows = []
    for item in sample:
        pid = item["id"]
        tier = item["tier"]
        if pid not in zh_reports or pid not in en_reports:
            continue
        z = extract_report(zh_reports[pid], "zh")
        e = extract_report(en_reports[pid], "en")
        for name, _, _ in ORGAN_SPECS:
            zu, eu = z[name], e[name]
            state_flag = decide_state_flag(zu, eu)
            over_norm = is_over_normalization(zu, eu)
            flag = "over_normalization" if over_norm else state_flag
            match = "Y" if zu["state"] == eu["state"] else "N"   # 四状态主标签一致性
            rows.append({
                "id": pid, "tier": tier, "organ": name,
                "zh_state": zu["state"], "en_state": eu["state"],
                "match": match,
                "zh_uptake_proposition": zu["proposition"],
                "en_uptake_proposition": eu["proposition"],
                "zh_unknown_reason": zu["unknown_reason"],
                "en_unknown_reason": eu["unknown_reason"],
                "trigger": (zu["trigger"] + " | " + eu["trigger"]).strip(" |"),
                "zh_trigger": zu["trigger"], "en_trigger": eu["trigger"],
                "evidence": "ZH: " + (zu["evidence"] or "—") + " || EN: " + (eu["evidence"] or "—"),
                "zh_evidence": zu["evidence"], "en_evidence": eu["evidence"],
                "zh_qualifiers": qualifiers_to_str(zu["qualifiers"]),
                "en_qualifiers": qualifiers_to_str(eu["qualifiers"]),
                "flag": flag, "state_flag": state_flag,
                "confidence": zu["confidence"] + " / " + eu["confidence"],
                "rule_version": RULE_VERSION,
            })
    return rows

ALIGN_COLS = ["id", "tier", "organ", "zh_state", "en_state", "match",
              "zh_uptake_proposition", "en_uptake_proposition",
              "zh_unknown_reason", "en_unknown_reason",
              "trigger", "zh_trigger", "en_trigger",
              "evidence", "zh_evidence", "en_evidence",
              "zh_qualifiers", "en_qualifiers",
              "flag", "state_flag", "confidence", "rule_version"]

def write_alignment(rows, path):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ALIGN_COLS)
        w.writeheader()
        w.writerows(rows)

def compute_stats(rows, sample, meta=None):
    meta = meta or {}
    n_pairs = len(set(r["id"] for r in rows))
    total_cells = len(rows) * 2
    zh_known = sum(r["zh_state"] != "Unknown" for r in rows)
    en_known = sum(r["en_state"] != "Unknown" for r in rows)
    known_cells = zh_known + en_known
    comparable = [r for r in rows if r["zh_state"] in KNOWN and r["en_state"] in KNOWN]
    agree = [r for r in comparable if r["zh_state"] == r["en_state"]]
    both_unknown = [r for r in rows if r["state_flag"] == "both_unknown"]

    by_organ = {}
    for name, _, _ in ORGAN_SPECS:
        rr = [r for r in rows if r["organ"] == name]
        comp = [r for r in rr if r["zh_state"] in KNOWN and r["en_state"] in KNOWN]
        ag = [r for r in comp if r["zh_state"] == r["en_state"]]
        by_organ[name] = {
            "cells": len(rr),
            "zh_known": sum(r["zh_state"] != "Unknown" for r in rr),
            "en_known": sum(r["en_state"] != "Unknown" for r in rr),
            "comparable": len(comp), "agree": len(ag),
            "agreement": round(len(ag) / len(comp), 3) if comp else None,
            "both_unknown": sum(r["state_flag"] == "both_unknown" for r in rr),
            "over_normalization": sum(r["flag"] == "over_normalization" for r in rr),
        }

    by_state = {}
    for st in ["Increased", "Decreased", "Absent", "Normal"]:
        zh_group = [r for r in comparable if r["zh_state"] == st]
        ag = [r for r in zh_group if r["zh_state"] == r["en_state"]]
        by_state[st] = {
            "zh_occurrences": sum(r["zh_state"] == st for r in rows),
            "en_occurrences": sum(r["en_state"] == st for r in rows),
            "zh_as_reference_comparable": len(zh_group),
            "zh_as_reference_agree": len(ag),
            "agreement": round(len(ag) / len(zh_group), 3) if zh_group else None,
        }

    flag_counts = defaultdict(int)
    state_flag_counts = defaultdict(int)
    for r in rows:
        flag_counts[r["flag"]] += 1
        state_flag_counts[r["state_flag"]] += 1

    # 修正2：四种 Unknown 成因分布（zh / en 分别统计，含已知状态空白计数）
    unknown_reason_dist = {}
    for lang in ("zh", "en"):
        c = defaultdict(int)
        for r in rows:
            reason = r[lang + "_unknown_reason"]
            c[reason if reason else "known(blank)"] += 1
        unknown_reason_dist[lang] = dict(c)

    # 修正1：摄取命题分布（zh / en）
    proposition_dist = {}
    for lang in ("zh", "en"):
        c = defaultdict(int)
        for r in rows:
            c[r[lang + "_uptake_proposition"]] += 1
        proposition_dist[lang] = dict(c)

    # 修正1：over_normalization 案例
    on_cases = [r for r in rows if r["flag"] == "over_normalization"]
    over_normalization = {
        "trigger_rule": "zh proposition in {not_increased, not_abnormal, physiological} AND en proposition == normal_explicit",
        "count": len(on_cases),
        "rate_of_all_cells": round(len(on_cases) / len(rows), 4) if rows else None,
        "rate_of_comparable": round(len(on_cases) / len(comparable), 4) if comparable else None,
        "cases": [{
            "id": r["id"], "tier": r["tier"], "organ": r["organ"],
            "zh_state": r["zh_state"], "en_state": r["en_state"],
            "zh_uptake_proposition": r["zh_uptake_proposition"],
            "en_uptake_proposition": r["en_uptake_proposition"],
            "zh_trigger": r["zh_trigger"], "en_trigger": r["en_trigger"],
            "zh_evidence": r["zh_evidence"], "en_evidence": r["en_evidence"],
        } for r in on_cases],
    }

    # 按癌种分层（meta_data.csv 的 Disease）
    by_disease = {}
    if meta:
        disease_of = {}
        for pid, dis in meta.items():
            disease_of[pid] = dis
        grouped = defaultdict(list)
        for r in rows:
            grouped[disease_of.get(r["id"], "UNKNOWN")].append(r)
        for dis, rr in sorted(grouped.items()):
            comp = [r for r in rr if r["zh_state"] in KNOWN and r["en_state"] in KNOWN]
            ag = [r for r in comp if r["zh_state"] == r["en_state"]]
            bu = [r for r in rr if r["state_flag"] == "both_unknown"]
            by_disease[dis] = {
                "n_pairs": len(set(r["id"] for r in rr)),
                "cells": len(rr), "comparable": len(comp), "agree": len(ag),
                "comparable_agreement": round(len(ag) / len(comp), 4) if comp else None,
                "including_both_unknown_rate": round((len(ag) + len(bu)) / len(rr), 4) if rr else None,
                "over_normalization": sum(r["flag"] == "over_normalization" for r in rr),
            }

    # 候选真实漂移（沿用预实验判读逻辑：两侧已知但不一致，且非规则伪影）
    artifact_hints = {
        "05": ["喉咽", "舌尖", "牙槽", "tongue", "laryngopharyn"],
        "06": ["牙槽", "喉咽", "舌尖", "alveolar", "tooth"],
        "08": ["淋巴结", "纵隔", "肺门", "腋下", "lymph node", "axillary", "mediastin"],
        "11": ["膈前组", "淋巴结", "diaphragm", "lymph node"],
        "14": ["胰尾", "实性部分", "pancrea"],
        "16": ["附件区", "adnexal"],
        "17": ["附件区", "淋巴结", "adnexal", "lymph node"],
        "19": ["鼻咽", "骨髓", "肩周", "淋巴结", "颅底", "nasopharyn", "skull base", "bone marrow"],
    }
    def assess(r):
        blob = r["zh_evidence"] + " " + r["en_evidence"]
        hits = [w for w in artifact_hints.get(r["organ"][:2], []) if w in blob]
        if hits:
            return "rule_artifact_suspected(跨器官词:" + "/".join(hits) + ")"
        return "candidate_real_drift(需人工复核)"

    candidate_drift = []
    for r in rows:
        if r["state_flag"] in ("negation_mismatch", "state_reversal", "uncertainty_shift",
                               "severity_shift", "state_mismatch"):
            a = assess(r)
            if a.startswith("candidate_real_drift"):
                candidate_drift.append({
                    "id": r["id"], "tier": r["tier"], "organ": r["organ"],
                    "zh_state": r["zh_state"], "en_state": r["en_state"],
                    "flag": r["state_flag"], "assessment": a,
                    "zh_evidence": r["zh_evidence"], "en_evidence": r["en_evidence"],
                })

    # 人工复核清单：所有 over_normalization + 所有 candidate_real_drift 的 id
    review_ids = sorted(set([c["id"] for c in on_cases] +
                            [c["id"] for c in candidate_drift]))

    stats = {
        "meta": {
            "generated_by": "organ_state_extractor.py v2.0 (rule-based, formal analysis)",
            "rule_version": RULE_VERSION,
            "ontology": "terminology_ontology_v01.md",
            "data_sources": ["petwb_data/report_zh_v3.csv", "petwb_data/report_en.csv",
                             "petwb_data/meta_data.csv" if meta else "p0_sample30.csv"],
            "n_pairs": n_pairs, "organs": len(ORGAN_SPECS),
            "states": ["Increased", "Decreased", "Absent", "Normal", "Unknown"],
            "propositions": PROPOSITIONS,
            "fixes_applied": ["修正1_uptake_proposition_over_normalization",
                              "修正2_四元unknown_reason", "修正3_五步分层与证据元信息",
                              "修正4_器官覆盖三层审计(仅full)"],
            "note": "四状态主分析；未提及保留为 Unknown，不自动改为 Normal。两个一致率口径固定且分开报告。",
        },
        "extraction_success": {
            "total_organ_language_cells": total_cells, "known_cells": known_cells,
            "overall_rate": round(known_cells / total_cells, 4) if total_cells else None,
            "zh_known": zh_known, "zh_rate": round(zh_known / len(rows), 4) if rows else None,
            "en_known": en_known, "en_rate": round(en_known / len(rows), 4) if rows else None,
            "both_unknown_cells": len(both_unknown),
        },
        "agreement": {
            "comparable_pairs": len(comparable), "agree_pairs": len(agree),
            "overall_comparable_rate": round(len(agree) / len(comparable), 4) if comparable else None,
            "including_both_unknown": {
                "numerator": len(agree) + len(both_unknown), "denominator": len(rows),
                "rate": round((len(agree) + len(both_unknown)) / len(rows), 4) if rows else None,
            },
            "by_organ": by_organ, "by_state": by_state,
        },
        "over_normalization": over_normalization,
        "unknown_reason_distribution": unknown_reason_dist,
        "proposition_distribution": proposition_dist,
        "flag_distribution": dict(flag_counts),
        "state_flag_distribution": dict(state_flag_counts),
        "candidate_real_drift": {"count": len(candidate_drift), "cases": candidate_drift},
        "manual_review": {
            "definition": "所有 over_normalization 与 candidate_real_drift 的 id",
            "n_ids": len(review_ids), "ids": review_ids,
        },
        "by_disease": by_disease,
        "limitations": [
            "规则按句级锚点定位器官，跨句/跨段的同一器官描述可能被拆散。",
            "器官级 Increased 优先：器官内病灶与'余…未见异常'共存时取 Increased，可能高估异常比例。",
            "Absent 仅识别'摄取缺失/无摄取'，器官术后缺如判为 Unknown。",
            "命题层为规则启发式，not_increased/not_abnormal 与 normal_explicit 的区分依赖固定表达模板，未覆盖所有译法。",
            "uncertain_in_report 仅在四状态主标签仍为 Unknown 时记录；若 hedge 描述同时命中状态词，则状态为已知、reason 留空（hedge 记录在限定词列）。",
            "英文报告为机器译文，存在 mojibake，但状态词未受影响。",
            "本结果为规则预标注，非人工金标准，不代表真实翻译错误率。",
        ],
    }
    return stats

# ---------------------------------------------------------------------------
# 八、修正4：器官覆盖三层审计
# ---------------------------------------------------------------------------
def build_audit_anchors():
    """① 原始词法锚点；② 本体同义词归一化后锚点。"""
    raw, syn = [], []
    for i, (_, zhs, ens) in enumerate(ORGAN_SPECS):
        raw.append(([re.compile(p, re.IGNORECASE) for p in zhs],
                    [re.compile(p, re.IGNORECASE) for p in ens]))
        zh_extra, en_extra = ONTOLOGY_SYNONYMS.get(i + 1, ([], []))
        zhs2 = list(zhs) + [re.escape(s) for s in zh_extra]
        ens2 = list(ens) + [re.escape(s) for s in en_extra]
        syn.append(([re.compile(p, re.IGNORECASE) for p in zhs2],
                    [re.compile(p, re.IGNORECASE) for p in ens2]))
    return raw, syn

def audit_coverage(zh_reports, en_reports):
    """修正4：对 490 对分别统计三层器官覆盖差异。

    ① 原始词法覆盖差异：原始锚点命中集合 zh vs en 不同的器官数
    ② 同义词归一化后覆盖差异：本体同义词扩展后仍不同的器官数
    ③ 正则标记需人工确认的差异：② 的差异中，锚点句含歧义/不确定标记者
    明确：只有 ③ 可能进入漂移结论；①② 仅反映词法/同义词覆盖，不等于漂移。
    """
    raw, syn = build_audit_anchors()
    amb_zh = [re.compile(p, re.IGNORECASE) for p in AMBIGUOUS_MARKERS["zh"]]
    amb_en = [re.compile(p, re.IGNORECASE) for p in AMBIGUOUS_MARKERS["en"]]
    org_names = [name for name, _, _ in ORGAN_SPECS]

    per_pair = []
    l1_by_organ = defaultdict(int)
    l2_by_organ = defaultdict(int)
    l3_by_organ = defaultdict(int)

    for pid in zh_reports:
        if pid not in en_reports:
            continue
        ztext, etext = normalize(zh_reports[pid]), normalize(en_reports[pid])
        zsent = split_sentences(ztext, "zh")
        esent = split_sentences(etext, "en")

        def present_sets(anchors):
            """anchors: 19 个 (zh_rx_list, en_rx_list)；返回命中的器官下标集合。"""
            zs = {i for i in range(len(org_names))
                  if any(rx.search(s) for rx in anchors[i][0] for s in zsent)}
            es = {i for i in range(len(org_names))
                  if any(rx.search(s) for rx in anchors[i][1] for s in esent)}
            return zs, es

        r_zs, r_es = present_sets(raw)
        s_zs, s_es = present_sets(syn)

        l1_diff = [org_names[i] for i in range(len(org_names)) if (i in r_zs) != (i in r_es)]
        l2_diff = [org_names[i] for i in range(len(org_names)) if (i in s_zs) != (i in s_es)]

        # ③ 对 ② 的差异，检查该器官锚点所在句是否含歧义/不确定标记
        l3_diff = []
        for i_name in l2_diff:
            i = org_names.index(i_name)
            z_hit = " ".join(s for s in zsent if any(rx.search(s) for rx in syn[i][0]))
            e_hit = " ".join(s for s in esent if any(rx.search(s) for rx in syn[i][1]))
            blob = z_hit + " " + e_hit
            if (any(rx.search(blob) for rx in amb_zh) or
                any(rx.search(blob) for rx in amb_en) or
                any(rx.search(blob) for rx in UNCERTAINTY_RX)):
                l3_diff.append(i_name)

        for o in l1_diff: l1_by_organ[o] += 1
        for o in l2_diff: l2_by_organ[o] += 1
        for o in l3_diff: l3_by_organ[o] += 1

        if l2_diff or l1_diff:
            per_pair.append({
                "id": pid,
                "L1_raw_diff": l1_diff, "L1_count": len(l1_diff),
                "L2_synonym_diff": l2_diff, "L2_count": len(l2_diff),
                "L3_manual_confirm": l3_diff, "L3_count": len(l3_diff),
            })

    return {
        "meta": {
            "description": "器官覆盖三层审计（修正4）：① 原始词法 ② 同义词归一化 ③ 正则标记需人工确认",
            "rule_version": RULE_VERSION,
            "n_pairs": len([p for p in zh_reports if p in en_reports]),
            "organs": len(org_names),
            "sources": {"zh": "petwb_data/report_zh_v3.csv", "en": "petwb_data/report_en.csv",
                        "synonyms": "terminology_ontology_v01.md 第一节"},
            "layers": {
                "L1_raw_lexical": "原始锚点命中集合 zh vs en 不同的器官数（词法层，不等于漂移）",
                "L2_synonym_normalized": "本体同义词扩展后仍不同的器官数（归一化层，不等于漂移）",
                "L3_manual_confirmation": "L2 差异中锚点句含歧义/不确定标记者；只有本层可能进入漂移结论",
            },
        },
        "summary": {
            "pairs_with_any_diff_L1": sum(1 for p in per_pair if p["L1_count"] > 0),
            "pairs_with_any_diff_L2": sum(1 for p in per_pair if p["L2_count"] > 0),
            "pairs_with_any_diff_L3": sum(1 for p in per_pair if p["L3_count"] > 0),
            "total_organ_diffs_L1": sum(p["L1_count"] for p in per_pair),
            "total_organ_diffs_L2": sum(p["L2_count"] for p in per_pair),
            "total_organ_diffs_L3": sum(p["L3_count"] for p in per_pair),
        },
        "by_organ": {
            o: {"L1_diff_pairs": l1_by_organ.get(o, 0),
                "L2_diff_pairs": l2_by_organ.get(o, 0),
                "L3_manual_confirm_pairs": l3_by_organ.get(o, 0)}
            for o in org_names
        },
        "per_pair": per_pair,
        "conclusion_scope": (
            "只有 L3（regex 标记需人工确认的差异）可能进入漂移结论；"
            "L1/L2 仅反映词法与同义词覆盖差异，不可直接判为翻译错误。"
            "所有 L3 案例仍需人工逐对核对后方可写入论文结果。"
        ),
        "limitations": [
            "词法覆盖为关键词命中，非语义等价判定；同义词表来自本体 v0.1，未穷尽真实译法。",
            "歧义/不确定标记为启发式正则，可能漏标或过标，L3 数量仅为 triage 规模估计。",
            "覆盖差异不等于状态漂移；本审计不产出漂移率，只界定需人工复核的范围。",
        ],
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true", help="对全部 490 对运行（正式分析）")
    args = ap.parse_args()

    zh_reports, en_reports = load_reports()
    meta = load_meta()
    if args.full:
        sample = [{"id": i, "tier": ""} for i in zh_reports if i in en_reports]
        align_path = os.path.join(BASE, "state_alignment490.csv")
        stats_path = os.path.join(BASE, "full_stats.json")
    else:
        sample = load_sample(os.path.join(BASE, "p0_sample30.csv"))
        align_path = os.path.join(BASE, "state_alignment30.csv")
        stats_path = os.path.join(BASE, "p0_pilot_stats.json")

    rows = run(sample, zh_reports, en_reports)
    write_alignment(rows, align_path)
    stats = compute_stats(rows, sample, meta=meta)
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)

    if args.full:
        audit = audit_coverage(zh_reports, en_reports)
        audit_path = os.path.join(BASE, "organ_coverage_audit.json")
        with open(audit_path, "w", encoding="utf-8") as f:
            json.dump(audit, f, ensure_ascii=False, indent=2)
        print("wrote:", audit_path)
        print("coverage_L3_pairs:", audit["summary"]["pairs_with_any_diff_L3"])

    print("pairs_run:", stats["meta"]["n_pairs"])
    print("alignment_rows:", len(rows))
    print("extraction_success_rate:", stats["extraction_success"]["overall_rate"],
          "zh:", stats["extraction_success"]["zh_rate"],
          "en:", stats["extraction_success"]["en_rate"])
    print("agreement_comparable:", stats["agreement"]["overall_comparable_rate"],
          "n=", stats["agreement"]["comparable_pairs"])
    print("agreement_incl_both_unknown:", stats["agreement"]["including_both_unknown"]["rate"])
    print("over_normalization:", stats["over_normalization"]["count"])
    print("flag_distribution:", stats["flag_distribution"])
    print("unknown_reasons:", stats["unknown_reason_distribution"])
    print("manual_review_ids:", stats["manual_review"]["n_ids"])
    print("wrote:", align_path)
    print("wrote:", stats_path)

if __name__ == "__main__":
    main()
