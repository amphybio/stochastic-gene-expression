#!/bin/bash
set -eu -o pipefail

montage -density 300 -geometry +0+0 fig_entropy_vs_fano_{H,I}.eps fig_entropy_vs_fano.eps
montage -density 300 -geometry +0+0 fig_{entropy,information}_vs_epsilon_1_mu_50.eps fig_entropy_vs_epsilon.eps
for n in {1..3}; do
    montage -density 300 -geometry +0+0 fig_entropy_vs_mu_${n}_{H,H_ON,I}.eps fig_distributions_${n}.eps fig_entropy_vs_mu_${n}.eps
done

zip - >entropy_figs.zip fig_entropy_vs_{fano,mu_{1..3},epsilon}.eps
