# shared/npu — NPU (Ascend) 兼容层

所有 NPU 适配代码集中在此，**按设备类型条件激活**：NPU 上应用全部 workaround，CUDA/CPU 上为 no-op 或返回高性能默认值，避免在 GPU 上损失性能（eager attention、关闭 MHA fastpath、关闭 nested tensor 等）。

## 使用方式

训练入口（`KPLM/tasks/main.py`）在 import torch 之后调用一次：

```python
from shared.npu import apply_global_patches
apply_global_patches()  # 自动检测设备，非 NPU 环境 no-op
```

模型类中选择 attention 实现（替代硬编码 `"eager"`）：

```python
from shared.npu import attn_implementation
config._attn_implementation = attn_implementation(self.device)  # NPU -> eager, CUDA -> sdpa
```

## 包含的 patch（仅 NPU 生效）

| patch | 原因 |
|-------|------|
| `suppress_torch_npu_warnings` | torch_npu 的路径属主警告无实际影响 |
| `disable_nested_tensor` | NPU 不支持 `torch._nested_tensor_from_mask` |
| `disable_mha_fastpath` | eval 下 MHA fastpath 在长序列上极慢（111s → 0.28s/batch） |
| `patch_torch_save` | NPU tensor 先转 CPU 再序列化，兼容 torch_npu 检查 |

另有 `use_additive_attn_mask()`：`finetune_model.py` 的 TransformerAdapter 在 NPU 上用 additive float mask 替代 `src_key_padding_mask`。
