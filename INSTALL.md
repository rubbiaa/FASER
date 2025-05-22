Here’s a complete INSTALL.md file with all steps combined into a single Markdown file:

⸻


# FASER Software Installation Guide

This guide provides step-by-step instructions to build and install the FASER software and its dependencies.

---

## 🧰 Prerequisites

You must have the following installed:

### Common Tools
- `git`
- `cmake`
- `make`
- C++17 compiler (`g++` or `clang++`)

### Libraries and Tools
- **Boost**
- **Autotools** (Automake, Autoconf, Libtool, M4, Perl)
- **ROOT 6.20+**

### Install on Ubuntu/Debian:

```bash
sudo apt update
sudo apt install build-essential git cmake automake autoconf libtool m4 perl \
                 libboost-all-dev

Install on macOS (with Homebrew):

brew install boost automake autoconf libtool cmake


⸻

📦 Clone the Repository

git clone https://github.com/rubbiaa/FASER.git
cd FASER


⸻

🔨 Build Instructions

All components are built with make from the root of the FASER directory.

⸻

1. Build CLHEP

make clhep

This will:
	•	Configure and compile CLHEP
	•	Install it in CLHEP-install/

⸻

2. Build RAVE

make rave

This will:
	•	Run ./configure with CLHEP and Boost
	•	Compile and install to rave-install/

On macOS, Boost paths from Homebrew (/opt/homebrew) are used automatically.

⸻

3. Source ROOT

Before building GenFit, you must source the ROOT environment:

source /path/to/ROOT/root_install/bin/thisroot.sh

Replace /path/to/ROOT with your actual ROOT installation path.

⸻

4. Build GenFit

make genfit

This will:
	•	Build GenFit
	•	Install it in GenFit-install/

⸻

🧹 Clean Build Environment

To completely clean all builds:

make clean

⸻

💬 Support

If you encounter problems, open an issue on GitHub:
👉 https://github.com/rubbiaa/FASER/issues

⸻


