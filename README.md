# Spatial Outbreak Detection
Evaluate spatio-temporal effect in predicting possible disease outbreaks.


## Project Information

4 steps:
- Generate simple data and implement MCMC using a spatially independent model
- Generate spatial data and implement MCMC using both spatially independent and dependent models 
- Evaluate the spatio-temporal effect in terms of MCMC estimation
- Implement MCMC on Maryland Covid data using spatial dependent model 

## Code structure
```
├── generate_data_simple.R
├── generate_data_spatial.R
├── ckpt
│   └── ccks
│       └── role
│           ├── best.pdparams   # 模型参数
│           ├── .ipynb_checkpoints
│           └── test_pred.json  # 预测结果
├── conf
│   └── ccks
│       ├── enum_tag.dict
│       ├── event_schema.json   # 事件元素对应的类型字典
│       ├── role_tag.dict
│       └── trigger_tag.dict
├── data
│   └── ccks
│       ├── dev.json
│       ├── enum
│       │   ├── dev.tsv
│       │   ├── test.tsv
│       │   └── train.tsv
│       ├── pre_submit
│       │   ├── test_handle.csv
│       │   └── test_row.csv
│       ├── raw
│       │   ├── ccks_task1_eval_data.txt    # 金融原始测试语料数据
│       │   ├── ccks_task1_train.txt    # 金融原始训练语料数据
│       │   ├── task1_train.txt    # 医疗诊断原始训练语料数据
│       │   └── task1_unlabeld_val.txt    # 医疗诊断原始测试语料数据
│       ├── role
│       │   ├── dev.tsv
│       │   ├── test.tsv
│       │   └── train.tsv
│       ├── sentence
│       │   ├── dev.json    # 逐句切分后的数据集
│       │   ├── test.json
│       │   └── train.json
│       ├── test.json   # ccks2021转成lic2021格式后的数据
│       ├── train.json
│       └── trigger
│           ├── dev.tsv
│           ├── test.tsv
│           └── train.tsv
├── data_prepare.py     # BIO标注
├── log
│   ├── endpoints.log
│   └── workerlog.0
├── post_process.py     # 处理成最终提交的格式:result.txt
├── __pycache__
│   └── utils.cpython-37.pyc
├── run_event_element_extraction.sh     ###### 总控脚本 ######
├── run_sequence_labeling.sh    # 序列标注训练和预测脚本
├── sequence_labeling.py    # 模型代码
├── submit
│   └── result.txt      ###### 最终的提交文件 ######
│   └── 医疗三元组.json      ###### 根据训练数集和预测结果所生成的医疗三元组 ######
└── utils.py
```

## Running environment
Please refer to```requirement.txt```file.

## Running Command



### Maryland Covid Data

**Model parameters**



## Evaluation

simple model: precision: 0.67928, recall: 0.70428, f1: 0.69156
spatial model: precision: 0.68166, recall: 0.81932, f1: 0.74418
