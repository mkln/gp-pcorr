#!/bin/bash
SESSION="sims"
tmux new-session -d -s $SESSION

N=200
K=15
block_size=$(( (N + K - 1) / K ))

for i in $(seq 0 $((K-1))); do
  s_start=$(( i * block_size + 1 ))
  s_end=$(( (i+1) * block_size ))
  if [ $s_end -gt $N ]; then s_end=$N; fi
  if [ $s_start -gt $N ]; then break; fi

  if [ $i -gt 0 ]; then
    tmux split-window -t $SESSION
    tmux select-layout -t $SESSION tiled
  fi

  tmux send-keys -t $SESSION "Rscript pmatern_reps_1_fitting.r $s_start $s_end" Enter
done

tmux attach -t $SESSION