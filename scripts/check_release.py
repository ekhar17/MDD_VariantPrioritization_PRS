#!/usr/bin/env python3
"""Validate syntax, source-only contents and portable paths without study data."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
files = sorted(p for p in root.rglob('*') if p.is_file() and '.git' not in p.relative_to(root).parts)
errors = []
scripts = [p for p in files if p.suffix in {'.R', '.sh', '.sbatch', '.py'}]
private = re.compile(r'/(?:Users|home)/[^/\s]+/|/(?:proj|vast/som|nas)/[^/\s]+/(?:users|home)/' + '|' + '/' + 'work/users/')
for path in files:
    rel = path.relative_to(root)
    if 'vendor' in rel.parts and path.name != 'SHA256SUMS':
        errors.append(f'Third-party or generated vendor file included: {rel}')
    if path.suffix in {'.rds','.bk','.bed','.bim','.fam','.bgen','.sample','.RData','.sif','.zip','.gz','.pdf','.png'}:
        errors.append(f'Data/binary artifact included: {rel}')
    if path.name == 'paths.local.sh':
        errors.append('Local configuration must be excluded from the upload.')
    if path.suffix not in {'.pyc'}:
        try: text = path.read_text()
        except UnicodeDecodeError:
            errors.append(f'Unexpected binary: {rel}'); continue
        if private.search(text): errors.append(f'Original personal/cluster path remains: {rel}')
        # Common credential formats; names of environment settings are allowed.
        if re.search(r'gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,}|AKIA[0-9A-Z]{16}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----',text):
            errors.append(f'Possible embedded credential: {rel}')
for path in scripts:
    if path.suffix in {'.sh','.sbatch'}:
        p = subprocess.run(['bash','-n',str(path)],capture_output=True,text=True)
        if p.returncode: errors.append(p.stderr.strip())
    elif path.suffix == '.py':
        try: compile(path.read_text(),str(path),'exec')
        except SyntaxError as e: errors.append(str(e))
r_files = [str(p.relative_to(root)) for p in scripts if p.suffix == '.R']
if not shutil.which('Rscript'):
    errors.append('Rscript is required for R syntax validation.')
else:
    expression = 'fs <- c(' + ','.join(json.dumps(p) for p in r_files) + '); for (f in fs) parse(f); cat(length(fs), "R files parsed\\n")'
    p = subprocess.run(['Rscript','--vanilla','-e',expression],cwd=root,capture_output=True,text=True)
    if p.returncode: errors.append(p.stderr.strip())
    else: print(p.stdout.strip())
manifest = json.loads((root/'docs/SCRIPT_PROVENANCE.json').read_text())
for item in manifest:
    path = root/item['path']
    if not path.is_file(): errors.append('Missing provenance target: '+item['path'])
    elif hashlib.sha256(path.read_bytes()).hexdigest() != item['release_sha256']:
        errors.append('Release hash changed without provenance update: '+item['path'])
for item in json.loads((root/'docs/RENAMING_MAP.json').read_text()):
    if not (root/item['current_path']).is_file():
        errors.append('Missing renamed file: '+item['current_path'])
print(f'{len(files)} source/documentation files; {len(scripts)} executable scripts; {len(manifest)} study-script provenance entries.')
if errors:
    print('\n'.join(errors),file=sys.stderr)
    sys.exit(1)
print('PASS: syntax, source-only contents, path scan, credential-pattern scan and provenance targets.')
