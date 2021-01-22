# Project
Tools for finding heterogeneous treatment effects (and means) based on partitioning the covariate/feature space via full cross-cuts and solved via greedy search. A typical usage would be analyzing and experiment to find the high-level subgroups (a coarse partition that is useful to humans) that differ in their estimated treatment effects. 

This package is inspired by, and uses ideas from, [Causal Tree](https://github.com/susanathey/causalTree) but aims to have the partition be more interpretable and have better accuracy. It is slower, though for high-level partitions this is usually not an issue.

This project is currently in an advanced prototype stage. Issues may still be found in common usage. Please create issues for these!

## Contributing

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.opensource.microsoft.com.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

## Trademarks

This project may contain trademarks or logos for projects, products, or services. Authorized use of Microsoft 
trademarks or logos is subject to and must follow 
[Microsoft's Trademark & Brand Guidelines](https://www.microsoft.com/en-us/legal/intellectualproperty/trademarks/usage/general).
Use of Microsoft trademarks or logos in modified versions of this project must not cause confusion or imply Microsoft sponsorship.
Any use of third-party trademarks or logos are subject to those third-party's policies.
