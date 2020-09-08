# Parallel HuMAnN v1.0 #
- Update HuMAnN3.0 with parallel computing
- Speed up to 2x-4x times
- compatible with HuMAnN 3.0 

## 说明 ##
- 基于HuMAnN 3.0版本优化，兼容HuMAnN 3.0版本
- 主要使用Python Multiprocessing优化大文件处理流程
- 优化部分文件IO过程
- 整体实现2-4倍的加速
- 移除HuMAnN 3中memory-use参数的支持，影响运行效率


## 如何使用？ ##
请参见[wiki](../../wiki)


## TODO ##
- [x] 继续优化Diamond Output post-processing  
- [x] 在更多数据集上测试结果一致性

## HUMAnN Src ##

----

 * For complete documentation, please see the [HUMAnN User Manual](http://huttenhower.sph.harvard.edu/humann/manual).

 * Please direct questions to the [HUMAnN google group](https://groups.google.com/forum/#!forum/humann-users) (subscribe to receive HUMAnN news).
----