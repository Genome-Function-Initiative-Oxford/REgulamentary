# MDV Tools: Installation & Usage Guide

## Installation

**MDVTools** is currently available only on [TestPyPI](https://test.pypi.org/).

Install with:

```bash
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ mdvtools
```

**Requirements:**  
You must have `htslib` installed, as `bgzip` and `tabix` are required dependencies.

---

## Creating a Visualization

The following Python code will create an MDV project from the output of the [regulamentary pipeline](#):

```python
from mdvtools.conversions import create_regulamentary_project_from_pipeline
from mdvtools.serverlite import serve_project

p = create_regulamentary_project_from_pipeline(
    # 1. Path for the MDV Project folder (created if it doesn't exist)
    "/path/to/myproject",
    # 2. Path to the regulamentary config YAML
    "/path/to/config.yaml",
    # 3. Path to the pipeline output
    "/path/to/output",
    # Path to the ATAC/DNase BigWig file
    atac_bw="input/bigwigs/astrocyte_DNase_chr21.bw"
)
p.set_editable(True)
serve_project(p)
```

This will generate a visualization and launch a local server on port **5050**.  
Open your browser and visit [http://localhost:5050](http://localhost:5050) to view it.

### Required Parameters

1. **Project folder:** Path to a directory for the MDV project (created if it doesn’t exist). You need write access to this folder.
2. **Config YAML:** The YAML config file used to run the regulamentary pipeline. This specifies locations of track files.
3. **Pipeline output:** Path to the pipeline’s output (should contain the `ATAC`, `merge`, and `CTCF` folders).

### Optional Parameters

- **atac_bw:** Path to the ATAC/DNase bigWig file used for peak calling. If the pipeline was run with bigWig files, the pipeline config value will be used and this parameter can be omitted.
- **peaks:** Specify which peaks to use (`merged`, `ATAC`, or `CTCF`). Default is `merged`.
- **genome:** Genome assembly (`hg38` or `mm10`). If `'remove_blacklist'` is in the config, the genome will be inferred. Default is `hg38`.
- **openchrom:** Method for open chromatin identification (`DNase` or `ATAC`). Used for labeling only. Default is `DNase`.

---

## Visualizing the Project

![mdv](images/mdv_regulamentary.jpg)

To visualize an existing MDV project:

```python
from mdvtools.mdvproject import MDVProject
from mdvtools.serverlite import serve_project

p = MDVProject("/path/to/myproject")
server_project(p)

```
By default, the local server runs on port 5050 (can be changed using the `port` parameter).  
Open your browser to [http://localhost:5050](http://localhost:5050).

You can also visualise a project at the command line with:-
```
python -m mdvtools.serverlite /path/to/project
```

**Default view:**  
- A table with regulatory element data  
- An interactive genome browser  
- A metaplot (deepTools heatmap)  

You can filter data via charts; the table and metaplot update accordingly.  
Click an element in the table or metaplot to zoom into it in the browser.

**Metaplot controls:**  
- Zoom: Mouse wheel  
- Pan: Hold middle/right mouse button and drag  
- Settings: Click the cog icon to adjust color, group, or sort elements

---

## Creating a Web Page

To generate a static version of your visualization suitable for web hosting (e.g., GitHub Pages):

```python
from mdvtools.mdvproject import MDVProject

p = MDVProject("/path/to/myproject")
p.convert_to_static_page("myproject")
```

A folder (`myproject`) will be created in your current directory.  
Move this folder to your web server’s public directory. It will be accessible at:

```
https://myserver.com/myproject
```

---

## Example Visualizations

Explore these public Regulamentary MDV visualizations:

| Cell Type                        | URL                                                                 |
|----------------------------------|---------------------------------------------------------------------|
| B-cells                          | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7586             |
| Cardiac muscle cell              | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7587             |
| Endothelial cell of umbilical vein | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7588           |
| Keratinocyte                     | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7589             |
| Natural killer cell              | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7590             |
| CD4-positive                     | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7591             |
| Fibroblast of lung               | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7592             |
| CD14-positive monocyte           | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7593             |
| Fibroblast of dermis             | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7594             |
| Osteoblast                       | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7595             |
| Skeletal muscle myoblast         | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7596             |
| Mammary epithelial cell          | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7597             |
| Astrocyte                        | https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7598             |

---