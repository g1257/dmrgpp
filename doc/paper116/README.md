## JOSS DMRG++

## Building

From the repository root, run:

```sh
pandoc dmrgpp.md --citeproc -o dmrgpp.pdf
```

Pandoc reads the bibliography specified by the `bibliography` field in
`dmrgpp.md` and uses `--citeproc` to resolve the citations and generate the
reference list.
