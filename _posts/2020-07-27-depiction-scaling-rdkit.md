---
layout: post
title: Depicting legacy structure drawings
subtitle: Scaling molecules with the RDKit 
tags: [rdkit, python, cheminformatics]
comments: true
---

I assume most reading this blog post have also read and applied the [blog posts](https://rdkit.blogspot.com/2015/02/new-drawing-code.html) from Greg Landrum about 
RDKits "not so new anymore" drawing code.  You may ask yourself why should I read another blog post about 2D depiction with the RDKit?  In this post I will focus on creating 2D depictions for example for display in web applications from real-world legacy systems. You know, system that have been feed with data for decades and went through multiple migrations possibly even between vendors and structure formats.

## Some Context

Or maybe you don't actually know that problem? Different industries, different systems. Let's just say not everyone has in-house tailored systems with tons of fancy structure normalization and correction mechanisms. Some rely on 3rd party vendors which have limited or no such mechanisms. What the user draws, is what is saved in the database. This not only includes orientation of the structures but more importantly it's style like label size, bond thickness and bond length.

