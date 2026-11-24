#!/usr/bin/env python3
"""Build a TEclass2 run config from the repo's base config + a model directory.

Usage:
    teclass2_make_config.py <base_config.yml> <input_fasta> <device> \
                            <model_dir> <out_config.yml> [--checkpoint NAME]

Takes the repo's config.yml as the architecture base, then:
  * enables the classify branch (classification: true),
  * points the input FASTA at our curated library (classification_dataset_path),
  * sets the compute device (cuda/cpu),
  * resolves the trained checkpoint inside <model_dir> and writes the three keys
    TEclass2 concatenates to load it:
        load_model(model_save_path + model_name + '/' + from_checkpoint)

The checkpoint is auto-discovered (a HuggingFace dir with config.json +
pytorch_model.bin/model.safetensors, else a *.bin/*.pt/*.pth/*.ckpt file). Pass
--checkpoint to override with a path relative to <model_dir> (or absolute) if the
archive layout differs.
"""
import argparse
import os
import sys
import yaml


def find_checkpoint(model_dir):
    """Return the path TEclass2 should load, or None."""
    hf_dirs, weight_files = [], []
    for root, _dirs, files in os.walk(model_dir):
        if 'config.json' in files and (
            'pytorch_model.bin' in files or 'model.safetensors' in files
        ):
            hf_dirs.append(root)
        for fname in files:
            if fname.endswith(('.bin', '.pt', '.pth', '.ckpt')):
                weight_files.append(os.path.join(root, fname))

    candidates = hf_dirs or weight_files
    if not candidates:
        return None

    # Prefer an explicit "best"/"final" checkpoint, else take the last sorted
    # (usually the highest epoch) for determinism.
    preferred = [c for c in candidates if any(
        k in os.path.basename(c).lower() for k in ('best', 'final')
    )]
    return sorted(preferred or candidates)[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('base_config')
    ap.add_argument('input_fasta')
    ap.add_argument('device')
    ap.add_argument('model_dir')
    ap.add_argument('out_config')
    ap.add_argument('--checkpoint', default=None,
                    help='checkpoint path relative to model_dir (or absolute)')
    args = ap.parse_args()

    with open(args.base_config) as fh:
        cfg = yaml.safe_load(fh) or {}

    cfg['classification'] = True
    cfg['classification_dataset_path'] = args.input_fasta
    cfg['te_db_path'] = args.input_fasta
    cfg['device'] = args.device

    if args.checkpoint:
        ckpt = args.checkpoint
        if not os.path.isabs(ckpt):
            ckpt = os.path.join(args.model_dir, ckpt)
    else:
        ckpt = find_checkpoint(args.model_dir)
        if not ckpt:
            sys.exit(
                "ERROR: no checkpoint found under {!r}. Inspect the extracted "
                "model and pass --checkpoint.".format(args.model_dir)
            )

    ckpt = os.path.abspath(ckpt)
    sys.stderr.write("TEclass2 checkpoint: {}\n".format(ckpt))

    # Reconstruct so that model_save_path + model_name + '/' + from_checkpoint == ckpt
    cfg['model_save_path'] = os.path.dirname(ckpt) + '/'
    cfg['model_name'] = ''
    cfg['from_checkpoint'] = os.path.basename(ckpt)

    with open(args.out_config, 'w') as fh:
        yaml.safe_dump(cfg, fh, default_flow_style=False)


if __name__ == '__main__':
    main()
