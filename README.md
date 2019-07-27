# Reproduce Fast Upsampling
> ### Reoproduced a paper:
> ### [*Fast image/video upsampling*][1]

![cover](https://github.com/zzeitt/Reproduce-FastUpsampling/blob/master/results/cover_comparison.png)
> Low-resolution image vs Upsampled image, with magnification factor being 4.

## Structure
* DEMO_color.m 👈 [entrance file]
* nbDeconv.m
* optMu.m

Other scripts include:

* DEMO_gray.m
* DEMO_gray_func.m
* pixeldup.m (Source: [https://searchcode.com/codesearch/view/29472841/][2] retrieved in July 2019.)
* fftCGSRaL.m & asetupLnormPrior.m (Source: [http://zoi.utia.cas.cz/deconv_sparsegrad][3] retrieved in July 2019.)

[1]: https://dl.acm.org/citation.cfm?id=1409106
[2]: https://searchcode.com/codesearch/view/29472841/
[3]: http://zoi.utia.cas.cz/deconv_sparsegrad