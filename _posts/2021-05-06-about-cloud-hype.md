---
layout: post
title: Rant about the cloud hype
subtitle: the right tool for the right job
tags: [cloud, hardware]
comments: true
---

I admittedly spend too much time on tech forums and doing so I have observed
everyone and their dog wanting to deploy their "stuff" on the cloud. 

> Hi there, I created this new awesome AI model with <insert cool new hyped tech>. How can I deploy it to the cloud?

If you don't already know how to, you probably should not because you lack the reason to do so. 

## It's about cost

I think the core issue here is that the general believe is that the cloud is cheap and local hardware is expensive. This simply isn't true for almost all use-cases and especially not compute-heavy workloads for which the cloud is extremely expensive and often doesn't actually offer an ideal configuration for your needs. 

#### Is hardware actually expensive?

Just go to a local shops website and configure a powerful system. I tried several. You can get a a 32-core machine with a powerful GPU (RTX 2080 Super), 128 GB of RAM and lots of SSD storage for less than $7k. For 10k you get 64-cores, 256 GB RAM and a RTX 2080 TI. You can adjust the cores, RAM and GPU power according to your needs but the point is one can get a huge amount of compute power for 10k or less. More than you will actually will be able to easily use without the right software.

If you are very tight on a budget, old, used server hardware can also be a great option as it usually is dirt cheap at the cost of being a bit slower and having higher power use. But that would be the budget option. 

The conclusion is that you can get a huge amount of compute power for relative small amounts of money compared to say salary costs or office rental costs.

#### Cost isn't only monetary

You looked at the Amazon prices and calculated according to projected usage the price for going with the cloud compared to buying hardware. You decide the cloud is cheaper and go with the cloud. Wait. What? What about usability? How do your data scientists get the data in the cloud? How to the update it? What's the policy when to update and when to use the cloud? Of course this assumes you aren't going all-in with dedicated servers on a VPC which with 100% guarantee cost you more than buying hardware.

Your calculations can simply be too optimistic. But even if the are about right it will make your core employees be hesitant to use your compute resources because each use has a cost attached. Is it actually worth it to try this? It can hinder creativity and new insights. It also means to log in, start your instance, wait till it is available and only then can you start transferring over new data which means more waiting. The other option is to open up your database so you can access it from the cloud. Now you add a huge security cost you forgot in your calculation. Either was without dedicated instances on a VPC, usability and/or security will suffer.

## What's your point?

The cloud is cool if you are operating a web application or web service that has large usage spikes depending on work hours, day of week or even season.  Covering the huge spikes with own hardware in these cases can be too costly and hence on-demand makes a lot of sense.

However if your scenario is "compute" and you can basically saturate as much hardware as you can get your hands on (AI training), then it's simply cheaper to buy the hardware yourself. If you are using it close to 24/7 at 100% the clouds on-demand advantage becomes irrelevant. 