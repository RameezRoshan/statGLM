import os
import sys
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether, HRFlowable
)
from reportlab.pdfgen import canvas

class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        super(NumberedCanvas, self).__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_header_footer(num_pages)
            super(NumberedCanvas, self).showPage()
        super(NumberedCanvas, self).save()

    def draw_header_footer(self, page_count):
        self.saveState()
        self.setFont("Helvetica", 8)
        self.setFillColor(colors.HexColor("#718096"))
        
        # Header (pages > 1)
        if self._pageNumber > 1:
            self.drawString(42, 750, "statGLM: Guide to Publishing on Git, GitHub, and CRAN")
            self.setStrokeColor(colors.HexColor("#E2E8F0"))
            self.setLineWidth(0.5)
            self.line(42, 744, letter[0] - 42, 744)

        # Footer (all pages)
        footer_text = f"Page {self._pageNumber} of {page_count}"
        self.drawRightString(letter[0] - 42, 36, footer_text)
        self.drawString(42, 36, "statGLM Package Deployment Documentation • Confidential & Professional")
        self.setStrokeColor(colors.HexColor("#E2E8F0"))
        self.setLineWidth(0.5)
        self.line(42, 46, letter[0] - 42, 46)

        self.restoreState()

def create_code_block(code_text, style_code):
    # Escape XML entities
    escaped = code_text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    formatted = "<br/>".join([line.replace(" ", "&nbsp;") for line in escaped.splitlines()])
    p = Paragraph(formatted, style_code)
    t = Table([[p]], colWidths=[letter[0] - 84])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor("#F7FAFC")),
        ('BOX', (0,0), (-1,-1), 0.75, colors.HexColor("#E2E8F0")),
        ('PADDING', (0,0), (-1,-1), 7),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ]))
    return t

def create_alert_box(title, text, style_body):
    content = f"<b><font color='#C53030'>[IMPORTANT] {title}</font></b><br/>{text}"
    p = Paragraph(content, style_body)
    t = Table([[p]], colWidths=[letter[0] - 84])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor("#FFF5F5")),
        ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#FEB2B2")),
        ('LINELEFT', (0,0), (0,0), 3.5, colors.HexColor("#E53E3E")),
        ('PADDING', (0,0), (-1,-1), 8),
    ]))
    return t

def build_pdf(filename="statGLM_Publishing_Guide.pdf"):
    doc = SimpleDocTemplate(
        filename,
        pagesize=letter,
        leftMargin=42,
        rightMargin=42,
        topMargin=42,
        bottomMargin=42
    )

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=20,
        leading=24,
        textColor=colors.HexColor("#1A365D"),
        spaceAfter=4
    )

    subtitle_style = ParagraphStyle(
        'DocSubTitle',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=11,
        leading=15,
        textColor=colors.HexColor("#4A5568"),
        spaceAfter=15
    )

    h1_style = ParagraphStyle(
        'H1',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=14,
        leading=18,
        textColor=colors.HexColor("#1A365D"),
        spaceBefore=14,
        spaceAfter=6,
        keepWithNext=True
    )

    h2_style = ParagraphStyle(
        'H2',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=11,
        leading=15,
        textColor=colors.HexColor("#2B6CB0"),
        spaceBefore=10,
        spaceAfter=4,
        keepWithNext=True
    )

    body_style = ParagraphStyle(
        'Body',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=9.5,
        leading=13.5,
        textColor=colors.HexColor("#2D3748"),
        spaceAfter=5
    )

    bullet_style = ParagraphStyle(
        'Bullet',
        parent=body_style,
        leftIndent=15,
        firstLineIndent=-10,
        spaceAfter=3
    )

    code_style = ParagraphStyle(
        'Code',
        parent=styles['Normal'],
        fontName='Courier',
        fontSize=8.5,
        leading=11.5,
        textColor=colors.HexColor("#1A202C")
    )

    story = []

    # Title block
    story.append(Paragraph("statGLM: Publication & Distribution Guide", title_style))
    story.append(Paragraph("Instructions for Publishing to Git/GitHub and Submitting to CRAN", subtitle_style))
    story.append(HRFlowable(width="100%", thickness=1.5, color=colors.HexColor("#2B6CB0"), spaceAfter=12))

    # Introduction summary badge
    intro_table = Table([[
        Paragraph("<b>Package:</b> statGLM | <b>Version:</b> 0.1.0 | <b>License:</b> MIT | <b>CRAN Status:</b> 0 errors ✔, 0 warnings ✔, 0 notes ✔", body_style)
    ]], colWidths=[letter[0] - 84])
    intro_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor("#EBF8FF")),
        ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BEE3F8")),
        ('PADDING', (0,0), (-1,-1), 6),
    ]))
    story.append(intro_table)
    story.append(Spacer(1, 4))

    # =========================================================
    # PART 1: Publishing on Git & GitHub
    # =========================================================
    story.append(Paragraph("Part 1: Publishing on Git &amp; GitHub", h1_style))
    story.append(Paragraph(
        "Publishing to GitHub enables public version tracking, community collaboration, issue tracking, and allows users worldwide to install the package directly via <code>remotes::install_github()</code>.",
        body_style
    ))
    story.append(Spacer(1, 4))

    # Step 1
    story.append(Paragraph("Step 1: Create an Empty Repository on GitHub", h2_style))
    story.append(Paragraph("1. Log into your GitHub account and navigate to <b>https://github.com/new</b>.", bullet_style))
    story.append(Paragraph("2. Set <b>Repository name</b> to: <code>statGLM</code>.", bullet_style))
    story.append(Paragraph("3. Set the repository visibility to <b>Public</b>.", bullet_style))
    story.append(Paragraph("4. <b>Important:</b> Leave 'Add a README file', 'Add .gitignore', and 'Choose a license' <b>unchecked</b> (the local repository already contains comprehensive, production-ready versions of all three files).", bullet_style))
    story.append(Paragraph("5. Click the green <b>Create repository</b> button.", bullet_style))
    story.append(Spacer(1, 4))

    # Step 2
    story.append(Paragraph("Step 2: Link the Local Repository &amp; Push", h2_style))
    story.append(Paragraph(
        "Open PowerShell or your terminal in the package root directory (<code>d:\\Antigravity\\Proc_GLM</code>) and run the following git commands:",
        body_style
    ))
    git_push_code = """# 1. Add your GitHub repository as the remote origin
git remote add origin https://github.com/<your-username>/statGLM.git

# 2. Ensure your default branch is named 'main'
git branch -M main

# 3. Push your commits and configure upstream tracking
git push -u origin main"""
    story.append(create_code_block(git_push_code, code_style))
    story.append(Paragraph("<i>Note: If you use SSH keys with GitHub, use <code>git@github.com:&lt;your-username&gt;/statGLM.git</code> instead.</i>", body_style))
    story.append(Spacer(1, 4))

    # Step 3
    story.append(Paragraph("Step 3: Create a Release / Tag on GitHub", h2_style))
    story.append(Paragraph(
        "Tagging releases marks milestones and lets users pin specific stable package versions (e.g. <code>v0.1.0</code>).",
        body_style
    ))
    tag_code = """# Create an annotated release tag locally and push to GitHub
git tag -a v0.1.0 -m "Release version 0.1.0: Initial release of statGLM"
git push origin v0.1.0"""
    story.append(create_code_block(tag_code, code_style))
    story.append(Paragraph(
        "Alternatively, navigate to your GitHub repository page, click <b>Releases</b> &gt; <b>Draft a new release</b>, specify tag <code>v0.1.0</code>, title it <code>statGLM 0.1.0</code>, paste the release summary from <code>NEWS.md</code>, and click <b>Publish release</b>.",
        body_style
    ))
    story.append(Spacer(1, 4))

    # Step 4
    story.append(Paragraph("Step 4: Enable Automated CI via GitHub Actions (Optional but Recommended)", h2_style))
    story.append(Paragraph(
        "You can configure GitHub to automatically run <code>R CMD check</code> on every commit across Ubuntu Linux, macOS, and Windows:",
        body_style
    ))
    ci_code = """# In R console:
usethis::use_github_action("check-standard")

# Commit and push the generated workflow:
# git add .github/workflows/R-CMD-check.yaml
# git commit -m "ci: add standard R-CMD-check GitHub Action workflow"
# git push"""
    story.append(create_code_block(ci_code, code_style))
    story.append(Spacer(1, 4))

    # User install from GitHub
    story.append(Paragraph("How Users Install the Package from GitHub", h2_style))
    story.append(Paragraph(
        "Once pushed to GitHub, users can install the package immediately using:",
        body_style
    ))
    user_gh_code = """# Install remotes if not yet installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install statGLM from GitHub
remotes::install_github("<your-username>/statGLM")"""
    story.append(create_code_block(user_gh_code, code_style))
    story.append(Spacer(1, 8))

    # =========================================================
    # PART 2: Publishing on CRAN
    # =========================================================
    story.append(Paragraph("Part 2: Submitting and Publishing on CRAN", h1_style))
    story.append(Paragraph(
        "CRAN (The Comprehensive R Archive Network) is the official R package repository. Packages on CRAN can be installed by any R user worldwide via the standard command <code>install.packages(\"statGLM\")</code>.",
        body_style
    ))
    story.append(Paragraph(
        "<b>Current Package Status:</b> The package has been verified under strict <code>R CMD check --as-cran</code> and achieves <b>0 errors, 0 warnings, 0 notes</b>.",
        body_style
    ))
    story.append(Spacer(1, 4))

    # Pre-submission checklist
    story.append(Paragraph("Pre-Submission Checklist", h2_style))
    story.append(Paragraph("1. <b>Maintainer Details in DESCRIPTION:</b> Update <code>Authors@R</code> with your actual name and an active, monitored email address where you can receive CRAN verification links.", bullet_style))
    story.append(Paragraph("2. <b>Submission Comments (cran-comments.md):</b> Already created and configured in <code>.Rbuildignore</code>, describing test environments and check status.", bullet_style))
    story.append(Paragraph("3. <b>Release Notes (NEWS.md):</b> Already created with the changelog for version 0.1.0.", bullet_style))
    story.append(Spacer(1, 4))

    story.append(create_alert_box(
        "Active Maintainer Email Required",
        "CRAN sends an automated confirmation email with a verification link immediately after submission. You MUST click this link within a few hours, otherwise your submission is automatically discarded.",
        body_style
    ))
    story.append(Spacer(1, 4))

    # Step 1: Remote Multi-Platform Verification
    story.append(Paragraph("Step 1: Remote Multi-Platform Verification", h2_style))
    story.append(Paragraph(
        "Before official submission, test the package against CRAN's dedicated Windows development servers and R-hub multi-OS builders:",
        body_style
    ))
    cran_check_code = """# In R console:
# Submit to CRAN Windows development server (win-builder)
devtools::check_win_devel()
devtools::check_win_release()

# (Optional) Check on R-hub Linux/macOS build servers
rhub::rhub_check()"""
    story.append(create_code_block(cran_check_code, code_style))
    story.append(Paragraph(
        "You will receive automated email reports from win-builder within 15–30 minutes confirming whether the package passed with zero errors, warnings, and notes.",
        body_style
    ))
    story.append(Spacer(1, 4))

    # Step 2: Build Production Source Tarball
    story.append(Paragraph("Step 2: Build the Production Source Tarball", h2_style))
    story.append(Paragraph(
        "Build the clean, compressed source distribution file (<code>.tar.gz</code>):",
        body_style
    ))
    build_code = """# In R console:
devtools::build()

# Output: statGLM_0.1.0.tar.gz created in the parent directory"""
    story.append(create_code_block(build_code, code_style))
    story.append(Spacer(1, 4))

    # Step 3: Submit to CRAN
    story.append(Paragraph("Step 3: Submit to CRAN", h2_style))
    story.append(Paragraph(
        "You can submit via either of the following two methods:",
        body_style
    ))
    story.append(Paragraph("<b>Method A: Automated via devtools (Recommended)</b>", body_style))
    submit_code = """# In R console:
devtools::submit_cran()"""
    story.append(create_code_block(submit_code, code_style))
    story.append(Spacer(1, 4))

    story.append(Paragraph("<b>Method B: Web Form Upload</b>", body_style))
    story.append(Paragraph("1. Open the CRAN submission portal in your browser: <b>https://cran.r-project.org/submit.html</b>", bullet_style))
    story.append(Paragraph("2. Enter your name and active email address.", bullet_style))
    story.append(Paragraph("3. Upload the built tarball: <code>statGLM_0.1.0.tar.gz</code>.", bullet_style))
    story.append(Paragraph("4. In the comments field, paste the contents of <code>cran-comments.md</code>.", bullet_style))
    story.append(Paragraph("5. Click <b>Upload Package</b>.", bullet_style))
    story.append(Spacer(1, 4))

    # Step 4: Confirm via Email
    story.append(Paragraph("Step 4: Confirm via Email (Crucial Action)", h2_style))
    story.append(Paragraph(
        "Within 5–15 minutes of submission, CRAN's submission system will send an email from <code>cran-submissions@r-project.org</code> with the subject <i>'CRAN submission statGLM 0.1.0'</i>.",
        body_style
    ))
    story.append(Paragraph("1. Open the email message immediately.", bullet_style))
    story.append(Paragraph("2. Click the verification URL provided inside the email.", bullet_style))
    story.append(Paragraph("3. A web page will load confirming: <i>'Submission of statGLM confirmed'</i>.", bullet_style))
    story.append(Spacer(1, 4))

    # Step 5: Review & Availability
    story.append(Paragraph("Step 5: Review &amp; Worldwide Availability", h2_style))
    story.append(Paragraph("After confirmation, the package undergoes CRAN's standard review workflow:", body_style))
    story.append(Paragraph("• <b>Automated Pre-Tests:</b> CRAN incoming servers run automated builds across Debian Linux, Fedora Linux, Windows, and macOS.", bullet_style))
    story.append(Paragraph("• <b>Manual Volunteer Inspection:</b> As a brand-new package (version 0.1.0), a CRAN volunteer team member inspects the description, license, and check logs.", bullet_style))
    story.append(Paragraph("• <b>Acceptance Notification:</b> Upon approval, you will receive an email stating: <i>'statGLM 0.1.0 is on its way to CRAN'</i>.", bullet_style))
    story.append(Paragraph("• <b>Mirror Distribution:</b> Within 24–48 hours, pre-compiled binary packages are generated and mirrored to worldwide CRAN servers.", bullet_style))
    story.append(Spacer(1, 4))

    # User install from CRAN
    story.append(Paragraph("How Users Install the Package from CRAN", h2_style))
    story.append(Paragraph(
        "Once accepted on CRAN, any user can install the package with standard R commands without needing devtools or Git:",
        body_style
    ))
    cran_install_code = """# Standard CRAN installation
install.packages("statGLM")

# Load and use
library(statGLM)"""
    story.append(create_code_block(cran_install_code, code_style))

    doc.build(story, canvasmaker=NumberedCanvas)
    print(f"Successfully generated {filename}")

if __name__ == '__main__':
    build_pdf()
