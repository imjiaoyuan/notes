# Gender Recognition

## 使用框架：

-   rembg
    - 去除叶片的背景
-   opencv
    - 将同一图片中的多个叶片分割开
    - 增加灰度处理，增强叶片特征
-   pillow
    - 将叶片处理成统一正方形尺寸的标准化数据
-   torchvision
    - 对标准化的图片进行数据增强
    - 添加黑色方块随即遮挡，防止过拟合

## 模型训练

-   Pytorch
    - 搭建、训练和优化神经网络模型
-   Torchvision
    - 预训练模型：**EfficientNet-B0**
        - 已学习丰富的通用视觉特征如纹理、边缘、形状等，只需在此基础上进行微调，使其适应小数据集

## 关键方法

迁移学习 (Transfer Learning)

```python
# 加载在 ImageNet 上预训练好的模型
model = models.efficientnet_b0(weights=models.EfficientNet_B0_Weights.IMAGENET1K_V1)
# 将模型的最终分类层替换为适应新任务的层
model.classifier[1] = nn.Linear(model.classifier[1].in_features, NUM_CLASSES)
```

学习率预热 (Learning Rate Warmup)

```python
if epoch < warmup_epochs:
    # 线性增加学习率
    lr_scale = (epoch + 1) / warmup_epochs
    for param_group in optimizer.param_groups:
        param_group['lr'] = LEARNING_RATE * lr_scale
```

余弦退火学习率 (Cosine Annealing LR)

```python
# 在 main 函数中定义调度器
scheduler = CosineAnnealingLR(optimizer, T_max=NUM_EPOCHS - WARMUP_EPOCHS, eta_min=1e-6)

# 在 train_one_epoch 函数的预热期后调用
if epoch >= warmup_epochs:
    scheduler.step()
```

AdamW 优化器

```python
optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE)
```

标签平滑 (Label Smoothing)

```python
criterion = nn.CrossEntropyLoss(label_smoothing=0.1)
```

模型检查点 (Model Checkpointing)

```python
if epoch_acc > best_acc:
    best_acc = epoch_acc
    torch.save(model.state_dict(), "best_model.pth")
```

## EfficientNet 模型

[EfficientNet: Rethinking Model Scaling for Convolutional Neural Networks](https://arxiv.org/abs/1905.11946)
