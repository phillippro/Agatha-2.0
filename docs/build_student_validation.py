"""Build the student validation protocol with Python and reportlab.

Companion training script: docs/examples/validation_smoke.R.
"""
from pathlib import Path
import math
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (
    BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, PageBreak,
    Table, TableStyle, Preformatted, Flowable, Image,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/pdf/agatha2_student_validation_protocol.pdf'
OUT.parent.mkdir(parents=True, exist_ok=True)
W, H = A4
M = 48
CW = W - 2*M
NAVY = colors.HexColor('#142F46')
TEAL = colors.HexColor('#007F86')
INK = colors.HexColor('#253847')
MUTED = colors.HexColor('#5C6C79')
LIGHT = colors.HexColor('#EDF4F6')
LINE = colors.HexColor('#D5E1E6')
ORANGE = colors.HexColor('#B96224')

styles = {
 'body': ParagraphStyle('body', fontName='Helvetica', fontSize=10.2,
   leading=14.3, textColor=INK, spaceAfter=8),
 'small': ParagraphStyle('small', fontName='Helvetica', fontSize=8.5,
   leading=11.6, textColor=MUTED, spaceAfter=6),
 'h1': ParagraphStyle('h1', fontName='Helvetica-Bold', fontSize=23,
   leading=27, textColor=NAVY, spaceAfter=15),
 'h2': ParagraphStyle('h2', fontName='Helvetica-Bold', fontSize=12,
   leading=16, textColor=TEAL, spaceBefore=9, spaceAfter=6, keepWithNext=True),
 'kicker': ParagraphStyle('kicker', fontName='Helvetica-Bold', fontSize=8.5,
   leading=12, textColor=TEAL, spaceAfter=8),
 'cell': ParagraphStyle('cell', fontName='Helvetica', fontSize=9,
   leading=12.2, textColor=INK),
 'th': ParagraphStyle('th', fontName='Helvetica-Bold', fontSize=9,
   leading=12, textColor=colors.white),
 'code': ParagraphStyle('code', fontName='Courier', fontSize=8.1,
   leading=10.6, textColor=INK, leftIndent=10, rightIndent=8,
   borderPadding=10, backColor=LIGHT, spaceBefore=4, spaceAfter=12),
}

story=[]
sections=[]

def p(text, style='body'):
    story.append(Paragraph(text, styles[style]))

def h(text): p(text,'h2')
def code(text): story.append(Preformatted(text.strip('\n'),styles['code']))
def gap(n=6): story.append(Spacer(1,n))

def table(headers, rows, widths=None):
    if headers[0]=='READING PATH':
        rows=[(label, f'<link href="#p{pages.split("-")[0]}" color="#007F86">{pages}</link>')
              for label,pages in rows]
    vals=[[Paragraph(x, styles['th']) for x in headers]]
    vals += [[Paragraph(str(x),styles['cell']) for x in row] for row in rows]
    t=Table(vals,colWidths=widths or [CW/len(headers)]*len(headers),hAlign='LEFT',repeatRows=1)
    t.setStyle(TableStyle([
      ('BACKGROUND',(0,0),(-1,0),NAVY),
      ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,LIGHT]),
      ('VALIGN',(0,0),(-1,-1),'TOP'),
      ('LEFTPADDING',(0,0),(-1,-1),9),('RIGHTPADDING',(0,0),(-1,-1),9),
      ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
      ('LINEBELOW',(0,0),(-1,0),0.5,NAVY),
      ('LINEBELOW',(0,-1),(-1,-1),0.5,LINE),
    ]))
    story.extend([t,Spacer(1,10)])

def note(title,text):
    t=Table([[Paragraph('<b>'+title+'</b><br/>'+text,styles['body'])]],colWidths=[CW])
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),LIGHT),
      ('BOX',(0,0),(-1,-1),0.5,LINE),('LINEBEFORE',(0,0),(0,0),3,TEAL),
      ('LEFTPADDING',(0,0),(-1,-1),12),('RIGHTPADDING',(0,0),(-1,-1),12),
      ('TOPPADDING',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
    story.extend([t,Spacer(1,10)])

def page(title, source):
    if story: story.append(PageBreak())
    n=len(sections)+1
    sections.append((title,source))
    p(f'<a name="p{n}"/>AGATHA2 / STUDENT VALIDATION PROTOCOL','kicker')
    p(title,'h1')

class CoverArt(Flowable):
    def __init__(self): super().__init__(); self.width=CW; self.height=170
    def draw(self):
        c=self.canv
        c.setFillColor(LIGHT); c.roundRect(0,0,CW,170,10,fill=1,stroke=0)
        c.setStrokeColor(LINE); c.setLineWidth(.6)
        for y in [35,85,135]: c.line(18,y,CW-18,y)
        c.setStrokeColor(TEAL); c.setLineWidth(2)
        path=c.beginPath()
        for i in range(480):
            x=18+i*(CW-36)/479
            y=85+37*math.sin(i/42)+11*math.sin(i/21+.3)
            if i==0:path.moveTo(x,y)
            else:path.lineTo(x,y)
        c.drawPath(path)
        for i in range(25):
            x=18+(i*19+7)%480*(CW-36)/479
            q=(x-18)*479/(CW-36)
            y=85+37*math.sin(q/42)+11*math.sin(q/21+.3)+5*math.sin(i*7)
            c.setStrokeColor(NAVY);c.setLineWidth(.6);c.line(x,y-7,x,y+7)
            c.setFillColor(NAVY);c.circle(x,y,2,stroke=0,fill=1)
        c.setFont('Helvetica',8);c.setFillColor(MUTED)
        c.drawString(18,13,'Illustrative signal and measurements - not a fitted data set')

class Workflow(Flowable):
    def __init__(self): super().__init__();self.width=CW;self.height=72
    def draw(self):
        c=self.canv
        labels=[('Prepare','validate tables'),('Search','signal + noise'),('Refine','orbit + MCMC'),('Diagnose','activity + time'),('Assess','stability')]
        bw=91; sp=(CW-5*bw)/4
        for i,(a,b) in enumerate(labels):
            x=i*(bw+sp)
            c.setFillColor(LIGHT);c.roundRect(x,10,bw,51,5,fill=1,stroke=0)
            c.setFillColor(NAVY);c.setFont('Helvetica-Bold',10);c.drawCentredString(x+bw/2,40,a)
            c.setFont('Helvetica',7.5);c.drawCentredString(x+bw/2,25,b)
            if i<4:
                c.setStrokeColor(TEAL);c.line(x+bw+2,35,x+bw+sp-2,35)

class GuideDoc(BaseDocTemplate):
    def afterFlowable(self, flow):
        if (isinstance(flow,Paragraph) and flow.style.name=='h1'
            and flow.getPlainText() in {title for title,_ in sections}):
            key=f'section-{self.page}'
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(flow.getPlainText(),key,0,False)

def footer(c,doc):
    n=doc.page
    c.saveState()
    if n>1:
        c.setFillColor(MUTED);c.setFont('Helvetica',8)
        c.drawString(M,H-32,'AGATHA2')
        c.drawRightString(W-M,H-32,'Protocol v1.0 | 03 October 2026 | Code 229dccc')
        c.setStrokeColor(LINE);c.line(M,H-39,W-M,H-39)
    c.setStrokeColor(LINE);c.line(M,45,W-M,45)
    c.setFillColor(MUTED);c.setFont('Helvetica',7)
    source=sections[n-1][1] if n<=len(sections) else 'Layout overflow'
    c.drawString(M,32,source)
    c.setFont('Helvetica-Bold',8);c.drawRightString(W-M,31,f'{n:02d}')
    c.restoreState()

# 1
page('Validate Agatha2', 'Student protocol v1.0 | Proposed experiments, not completed validation results')
gap(15)
p('Student experimental protocol','h1')
p('Detection reliability, orbital recovery and uncertainty calibration','h2')
p('This report specifies a reproducible test campaign for Agatha2. It combines observed benchmark stars, controlled simulations and injections into realistic time series. The deliverable is an evidence-based assessment of where the software works, where it fails and which conclusions remain conditional.')
gap(10); story.append(Workflow()); gap(14)
table(['DOCUMENT CONTROL','VALUE'],[
 ('Prepared','03 October 2026'),
 ('Code to reproduce','Branch multiset-fitting; commit 229dccc33dd0107ca765be8cd97a30c4aeb3f684'),
 ('Scope','BFP detection, circular/Keplerian fits, noise sensitivity, multi-set behavior, MCMC and diagnosis. Short N-body checks are an extension.'),
 ('Core sample','Four observed hosts; 15 challenge-file candidates; 1000 null and 2000 one-planet synthetic realizations.'),
 ('Starting exercise','Rscript docs/examples/validation_smoke.R'),
], [121,CW-121])
note('What counts as success','A complete, reproducible campaign is successful even when it exposes a software limitation. Do not alter data or thresholds to obtain a desired planet count. Separate detection, local parameter refinement and interpretation of stellar activity.')
p('The code revision is a local snapshot and was previously unpushed. The project lead must distribute that revision and these companion files to all students. A fresh GitHub clone may not contain them. No full validation campaign has been run for this report.','small')

# 2
page('Assignment and milestones','Protocol design | Roles can be combined for a single-student project')
table(['STAGE','STUDENT TASK','EXIT PRODUCT'],[
 ('A / Setup','Record code and environment; run the training exercise and relevant existing tests.','Environment record, logs and two training results.'),
 ('B / Audit','Inspect benchmark inputs and challenge provenance. Define units, cuts and instrument segments.','Input manifest and signed-off cleaning rules.'),
 ('C / Pilot','Run 100 null and 200 injected cases using pilot-only seeds. Measure runtime and check scoring.','Pilot report; frozen settings and compute estimate.'),
 ('D / Frozen test','Run independent seeds for 1000 null and 2000 injected cases. Analyze observed benchmarks separately.','All result rows, including failures; no truth-driven retuning.'),
 ('E / Calibration','Run preselected MCMC cases and stress experiments. Score once outputs are frozen.','Recovery, false-detection, interval and runtime summaries.'),
 ('F / Review','A second reader reproduces selected successes and failures from the saved inputs.','Final report, figures, machine-readable tables and rerun commands.'),
], [77,219,CW-296])
h('Suggested responsibilities')
p('<b>Project lead:</b> fixes the scope, supplies the exact revision, approves the pilot settings and retains hidden truth. <b>Data student:</b> audits inputs and generates cases. <b>Analysis student:</b> runs the frozen pipeline without accessing validation truth. <b>Reviewer:</b> checks matching, denominators, failed jobs and selected reruns. With one student, freeze the analysis outputs before opening truth files.')
h('Suggested schedule')
p('Week 1: setup, input audit and pilot. Week 2: observed benchmarks and the main deterministic campaign. Week 3: MCMC calibration and targeted stress cases. Week 4: independent reruns and final report. These are planning milestones; the pilot timing determines the actual schedule, especially for dense GP fits.')
table(['READING PATH','PAGES'],[
 ('Environment, targets and input rules','3-5'),
 ('Simulation design, campaign and runner','6-8'),
 ('Observed targets and scoring','9-10'),
 ('MCMC and stress experiments','11-12'),
 ('Records, deliverables and completion criteria','13-14'),
 ('Troubleshooting and reference material','15-16'),
], [CW-80,80])

# 3
page('Set up the environment','Sources: server.R; mcmc_func.R; periodoframe.R; docs/examples/validation_smoke.R')
p('Use the project lead\'s checkout and run commands from its root. The current source imports doMC and targets a compatible macOS/Linux environment. Record the R version and all package versions; do not assume that upgrading R or a package leaves numerical results unchanged.')
code('''git rev-parse HEAD
git status --short
Rscript -e 'sessionInfo()'
Rscript docs/examples/validation_smoke.R''')
p('The training script creates a new timestamped directory under results/student_validation/. It simulates 80 observations with a known 17.3-day circular signal, alongside a paired no-planet series. It saves inputs, settings, hashes, comparisons, figures, warnings and session information. This is a setup check; one null realization cannot calibrate false detections.')
h('Install only missing R dependencies')
code('''pkgs <- c("shiny", "minpack.lm", "magicaxis", "foreach",
          "doMC", "ramify", "fields", "MASS", "Matrix",
          "doParallel", "e1071", "kernlab", "glasso",
          "JPEN", "matrixcalc", "mvtnorm", "rootSolve")
missing <- pkgs[!vapply(pkgs, requireNamespace,
                       logical(1), quietly=TRUE)]
if (length(missing))
  install.packages(missing, repos="https://cloud.r-project.org")''')
h('Run the existing regression tests')
code('''Rscript tests/test_eprior.R
Rscript tests/test_fourier_kepler.R
Rscript tests/test_sequential_residual.R''')
p('Then run test_multiset_periodogram.R, test_multiset_mcmc.R, test_gp_periodogram.R, test_diagnosis.R, test_moving_periodogram.R and test_ptmcmc.R for the features in your campaign. Capture stdout, stderr, return codes and elapsed times. Stop to investigate failing tests before interpreting new experiments; do not edit tests simply to make them pass.')
note('What the existing tests establish','These checks cover selected numerical and functional behaviors, mostly on fixed seeds. They do not measure population-level false-positive rates, completeness or interval calibration. Those are the purpose of this protocol.')

# 4
page('Observed benchmark sample','Local file audit | Literature references [1-6] are listed on page 16')
p('The following are recommended benchmark roles. Published periods are comparison labels, not inputs to the blind period search. A selected archive subset may not support every published detection.')
table(['TARGET','SCIENTIFIC ROLE','LOCAL INPUT / LIMIT'],[
 ('GJ 436','Strong short-period signal near 2.64 days [1]; detection and orbit refinement.','GJ436_DACE_HARPS03.dat: 171 rows, 1532-day span.'),
 ('HD 40307','Established inner periods near 4.3, 9.6 and 20.4 days [2]; sequential detection.','HD40307_DACE_HARPS03.dat: 451 rows, 4138-day span.'),
 ('HD 189567','Signals near 14.3 and 33.7 days [3]; noise, activity and alias sensitivity.','HD189567_DACE_HARPS03.dat: 632 rows, 4209-day span.'),
 ('HD 147379','Multi-instrument consistency and red noise [4]. Analyze separately and jointly.','HD147379_{carmenes,keck,sophie}.dat: 114, 30 and 113 rows; spans 602, 5166 and 3358 days.'),
 ('HD 156846','Optional high-eccentricity case: P about 359.5 days and e about 0.85 [5].','Local HARPS_TERRA file has 112 rows but only a 14-day span. Obtain a full-orbit data set before using it for recovery.'),
 ('HD 166435','Activity-origin variability near 3.8 days [6]; test diagnostic flags.','Acquire the published RVs and activity data. This is an activity-origin signal benchmark, not proof that the star has no planets at any period.'),
], [81,199,CW-280])
h('Fitting challenge and solar extension')
p('The repository has 15 ChallengeDataSet*_HARPS.dat files. Match them to the original RV Fitting Challenge definitions and truth tables before scoring [7]. Published solutions make these reproducibility benchmarks; use newly generated hidden cases for an independent blind test.')
p('Optional: HARPS-N solar RVs and indicators [8] provide realistic stellar and instrumental variability. Audit the applied Solar System corrections and the time sampling before using them as a noise baseline. Preserve the original release and record any resampling.')
note('Use the actual data version','Counts above are rows in the local files, not numbers of independent nights or validated measurements. The three DACE files were fetched on 7 May 2026; their times are DACE rjd and their RVs were median-centered per target. Keep that provenance with every result.')

# 5
page('Audit and freeze the inputs','Sources: data/DACE_TEST_DATA.md; functions.R:read.agatha.table')
p('For each file, create one manifest row before analysis. Keep an untouched original and a cleaned copy; every removed row must have a recorded reason independent of the desired periodogram result.')
table(['RECORD','REQUIRED CHECK'],[
 ('Identity','Target, instrument/segment, source URL or local provenance, release date, filename and checksum.'),
 ('Units','Time system and zero point; RV/error units; activity-index definitions. Convert explicitly and record the conversion.'),
 ('Sampling','Raw and retained rows; unique nights; time span; gaps; median uncertainty; observations per season.'),
 ('Validity','Finite times/RVs, strictly positive finite errors, duplicated timestamps, constant columns, missing-value sentinels and ordering.'),
 ('Segmentation','Known hardware or reduction changes. Predefine offsets/segments from metadata, not from the best-fitting planet solution.'),
 ('Activity','Retain valid indicators for diagnosis. Do not turn an undefined indicator into a random series or silently interpret -1 as a real measurement.'),
], [93,CW-93])
code('''tab <- read.table("data/GJ436_DACE_HARPS03.dat", header=TRUE)
tab <- tab[order(tab[,1]), ]
stopifnot(all(is.finite(as.matrix(tab[,1:3]))),
          all(tab[,3]>0))
sum(duplicated(tab[,1]))
c(N=nrow(tab), span=diff(range(tab[,1])),
  median_error=median(tab[,3]))
tools::md5sum("data/GJ436_DACE_HARPS03.dat")''')
note('Bypass surprising loader transformations','read.agatha.table() sorts and drops repeated timestamps, inserts unit errors for a two-column table, and replaces unusable extra columns with random values. For validation, load and validate tables explicitly, supply the required three columns, and preprocess each proxy using its documented missing-data rule.')
h('Freeze the comparison before seeing the answer')
p('Create a benchmark reference table: required candidate periods, the specific published data/solution being compared, and whether the test is a blind search or local refinement. Do not require the same amplitude or uncertainty when using a different instrument subset, noise model, prior or time coverage.')
p('Save the generation seed, input hash and settings hash for every simulated case. Distinguish three partitions: training, pilot/tuning and untouched validation. Do not reuse pilot noise realizations in the final campaign.')

# 6
page('Generate controlled time series','Proposed generators | Validate simulation code independently before use')
p('Use two cadence templates: the validated GJ 436 and HD 40307 timestamps and quoted errors. Subtract each template\'s first time for simulation, retaining the original epoch in metadata. Let T be its span and s the median quoted error. Use the same generated realization for paired comparisons between fitting methods.')
code('''y_i = gamma_set(i) + beta * (t_i - t_ref)
      + sum_j RV_j(t_i) + r_i + epsilon_i
epsilon_i ~ Normal(0, dy_i^2 + jitter^2)
RV_j(t) = K_j * [cos(nu_j(t) + omega_j)
                 + e_j*cos(omega_j)]''')
p('Compute the true anomaly nu from Kepler\'s equation using an independently checked solver. Draw phase and argument of periastron uniformly on [0, 2*pi); for e=0, score phase rather than the individually degenerate angles. Set jitter=0.5s; use gamma=0 and beta=0 in the core campaign. Offsets and trends enter the stress tests.')
table(['NOISE CLASS','GENERATOR FOR THE CORE CAMPAIGN'],[
 ('W','r=0; retain heteroscedastic quoted errors and the injected white jitter.'),
 ('MA-like','Independent z_i~N(0,1). r_i=A*(z_i+theta_i*z_(i-1))/sqrt(1+theta_i^2), theta_i=0.6*exp(-dt_i/10 d), r_1=A*z_1.'),
 ('AR-like','r_1=A*z_1; r_i=rho_i*r_(i-1)+A*sqrt(1-rho_i^2)*z_i; rho_i=0.8*exp(-dt_i/10 d).'),
 ('SHO GP','Draw a zero-mean GP with oscillation period 30 d and damping time 60 d. Use Q=pi*tau/Prot and normalize covariance diagonal to A^2.'),
 ('Quasi-periodic GP','k(dt)=A^2*exp[-dt^2/(2*tau^2)-2*sin^2(pi*dt/Prot)/ell^2]; Prot=30 d, tau=60 d, ell=1. This is a deliberate mismatch to the fitted SHO kernel.'),
], [90,CW-90])
p('Use A=2s for the four correlated classes. The lag generators above are explicitly defined benchmarks; do not describe them as exact matches to every Agatha2 AR/MA implementation. For GP draws, use a positive-definite covariance and record any numerical diagonal regularization.','small')
note('Generator acceptance check','Verify zero-planet data contain no injected RV term. On noiseless tests verify period, amplitude and eccentricity against an independent orbit calculation. Check covariance symmetry, diagonal normalization and sample statistics over many draws; do not validate a likelihood only against its own generator.')

# 7
page('The frozen experiment matrix','Proposed allocation | Realizations are shared across methods, not counted repeatedly')
table(['CAMPAIGN','ALLOCATION','TOTAL'],[
 ('Pilot null','2 cadences x 5 noise classes x 10 seeds','100'),
 ('Pilot injection','2 cadences x 5 noise classes x 5 periods x 4 amplitudes x 1 seed','200'),
 ('Validation null','2 cadences x 5 noise classes x 100 fresh seeds','1000'),
 ('Validation injection','2 cadences x 5 noise classes x 5 periods x 4 amplitudes x 10 fresh seeds','2000'),
], [107,319,CW-426])
h('Parameter grid')
p('Use P/T = 0.01, 0.05, 0.2, 0.7 and 1.2, and K/s = 0.5, 1, 2 and 4. Randomize phase and noise independently for each realization. Cycle eccentricity through 0, 0.3, 0.6 and 0.85, rotating the assignment across cells so the 2000 injections contain 500 of each eccentricity. This is a balanced design overall, not a full factorial eccentricity grid.')
p('The 2000 injections are 200 cadence/noise/period/amplitude cells with 10 repeats each. Fine-cell rates will be noisy: publish binomial intervals and pool only according to rules fixed in advance. Cases with P/T=1.2 test partial-orbit behavior; report them separately from well-covered recovery.')
h('Analysis branches and cost')
p('Run every core realization with W, MA(1) and AR(1) fits, each with Nh=2, Keplerian refinement and MCMC off. Designate MA(1) as the primary pipeline for this proposed protocol; W and AR(1) are paired comparisons. This gives 9000 deterministic analyses, not 9000 independent realizations.')
p('For the GP comparison, preselect null repeats 1-20 in every cadence/noise cell (200 cases), and injection repeats 1-2 in every cell (400 cases). Fit these 600 cases with a free SHO GP as an additional branch. If the pilot shows the budget is excessive, reduce a balanced design before opening validation truth and record the revised allocation.')
note('Freeze before launch','Freeze the input manifest, seeds, selected cases, model branches, search range, ofac, thresholds, failure policy and matching rule. Changes motivated by validation results require a new version and fresh holdout seeds. Never select the best noise model after checking injected truth.')
p('Estimate total CPU time from measured pilot times by branch and data size, including long tails and failures. Limit concurrency to allocated resources, explicitly set Ncores, and cap each job using a timeout selected from the pilot. A timeout is a recorded outcome, not a row to delete.','small')

# 8
page('Build the batch analysis','Sources: functions.R:calc.1Dper; companion validation_smoke.R')
p('Use the R API rather than the legacy agatha2.R CLI. Start from the companion training runner, then replace its generated table and narrow training grid with your case loader and frozen configuration. The example below analyzes one already validated three-column table called tab.')
code('''suppressPackageStartupMessages({
  source("periodograms.R")
  source("periodoframe.R")
  source("functions.R")
})
renew <- TRUE; Nsamp <- 1
Tspan <- diff(range(tab[,1]))
data <- list(case=tab)
pp <- list(ns="RV", ofac=2, frange=c(1/(2*Tspan),1/2),
  per.type="BFP", per.target="case", sequence=TRUE,
  Nsig.max=3, per.type.seq="BFP", Niter=0,
  SigType="kepler", Nh=2, Nmas=1, Nars=0,
  Inds=list(0), lnBF.min=5, noise.model="ARMA")
out <- calc.1Dper(1, "RV", pp, data, Ncores=1)
out$model.comp$RV
out$Nopt$RV''')
table(['BRANCH','CHANGE FROM THE PRIMARY CONFIGURATION'],[
 ('White','Nmas=0; Nars=0; noise.model="ARMA".'),
 ('MA(1), primary','Nmas=1; Nars=0; noise.model="ARMA".'),
 ('AR(1)','Nmas=0; Nars=1; noise.model="ARMA".'),
 ('Free GP, subset','noise.model="GP"; gp.Prot=NA; gp.tau=NA; Niter=0. Record warnings and fitted hyperparameters.'),
], [114,CW-114])
p('The common search spans 2 days to 2T; frange is in frequency units. Do not narrow it around an injected or published period for the blind search. Use a denser grid on predeclared pilot cases to check that ofac=2 is adequate before freezing it.')
note('Each job must leave an audit trail','Save the complete out object, pp, input hash, random seed, start/end time, warnings, process exit status and model.comp. Catch errors and write a failure row. A plotting failure should not erase a completed numerical fit; record analysis and export status separately in the full runner.')
p('Training code is included and tested. The complete campaign scheduler, generators and scoring program remain student implementation tasks, to be verified during the pilot.','small')

# 9
page('Analyze observed benchmarks','Observed-data procedure | Do not treat a literature fit as an exact injected truth')
table(['STEP','REQUIRED ACTION'],[
 ('1 / Lock reference','Select a literature solution and its data version. Record expected periods and the scope of the available local subset before analysis.'),
 ('2 / Inspect','Plot RV/error against time, activity indicators, the sampling window and seasonal coverage. Document outliers using frozen quality rules.'),
 ('3 / Blind search','Apply the common 2-day to 2T search, W/MA/AR branches, Nh=2 and the sequential threshold. Do not supply the reference period as an optimizer start for this stage.'),
 ('4 / Compare','Match candidates to reference periods; list missing reference signals, extra signals, alias alternatives and the accepted count. Explain data-set differences.'),
 ('5 / Refine','Refit plausible candidates with MCMC where supported. Label any analysis initialized at a reference period as local refinement, not blind discovery.'),
 ('6 / Diagnose','Compare activity power, individual instruments and moving windows. Record the exact diagnosis rules and show plots for important disagreements.'),
 ('7 / Export','Save model comparisons, posterior summaries, residuals, figures, settings and logs. Write a short conclusion for every target, including failures.'),
], [82,CW-82])
h('HD 147379: the multi-set assignment')
p('Analyze each instrument separately, then all three together, then leave one instrument out at a time. Keep a common time system and consistent RV units, with one fitted offset per set. For each run report period, amplitude, uncertainty, offsets, trend and the accepted signal count. A sparse instrument need not independently detect a real combined signal.')
p('Add a controlled test: use these three actual time/error arrays and generate one shared orbit, known offsets and per-set noise. This supplies exact truth for a capability that observed data alone cannot calibrate. Include overlapping and non-overlapping epoch subsets to expose offset/trend degeneracy.')
h('Activity and real-noise injections')
p('For HD 166435, assess whether activity and time-consistency diagnostics flag the known activity-origin periodicity; a BFP peak by itself is not a failure. For solar data or residual injections, retain the original variability, add a known orbit and compare with the un-injected baseline using identical settings. Report these as robustness experiments, separate from exact-null false-positive calibration.')
note('No post-hoc cleaning','If a recovered peak changes after removing data, show both outcomes and the reason for the cut. Do not delete a season or activity index because it makes the literature planet easier to recover.')

# 10
page('Score detections and failures','Proposed scoring rules | Freeze tolerances before evaluating validation truth')
h('Candidate selection and matching')
p('Use model.comp$RV and Nopt$RV, not the existence of a phase plot, to determine acceptance. In this snapshot, the first fitted signal can remain in parameter products even when its comparison row is rejected. A missing/nonfinite lnBF is an invalid significance result for this campaign; do not inherit the generic code\'s automatic acceptance of an NA comparison.')
code('''f_true = 1 / P_true;  f_hat = 1 / P_hat
match if abs(f_hat - f_true) <= min(0.1*f_true, 0.5/T)''')
p('This operational tolerance requires agreement both at 10% in frequency and within half a baseline frequency resolution. Match accepted candidates to injections one-to-one, maximizing the number matched and then minimizing total frequency error. Do not count P/2, 2P, daily or annual aliases as true recovery; label these separately. For P&gt;T also report prediction error and broad period intervals.')
table(['METRIC','DEFINITION / DENOMINATOR'],[
 ('Completeness','Matched injected planets / all injected planets. A failed analysis yields no recovery; also show the rate conditional on completed analyses.'),
 ('Null false-positive rate','Null realizations with at least one accepted candidate / completed null realizations. Always report failed fraction and all-attempt bounds: F/M to (F+U)/M, with F false-positive jobs, U failed jobs, M attempted jobs.'),
 ('Extra detections','Unmatched accepted candidates per completed realization; also report their fraction among accepted candidates.'),
 ('Multiplicity','Confusion table of true planet count versus accepted count. Add a distinct failed-job category.'),
 ('Parameter errors','Median and spread of P_hat-P_true, K_hat-K_true and e_hat-e_true for correctly matched signals. Keep undetected and mismatched cases visible in separate totals.'),
 ('Runtime / reliability','Median and upper-tail elapsed time, memory if measured, warnings, errors and timeouts by branch and data size.'),
], [111,CW-111])
p('Give 95% binomial intervals for rates. For pooled rates over heterogeneous cells, describe the design weights; use resampling by whole realization for stratified summaries and paired method differences. Three noise fits to one realization are not three independent observations.','small')
note('Threshold calibration','Keep ln(BF)&gt;5 as the initial operational threshold; it is not a known false-alarm probability. A proposed goal is a null false-positive rate below 1%, assessed with uncertainty. Zero events among 1000 completed comparable null cases gives only an approximately 0.3% one-sided 95% upper bound, not a zero true rate.')

# 11
page('Calibrate MCMC uncertainties','Sources: functions.R:mcfit, mcfit.multiset; mcmc_func.R:run.ptmcmc; references [9,10]')
p('Preselect 200 injection cases before seeing results: 20 from each cadence/noise combination, stratified over amplitude, period and eccentricity. Run the supported single-set ARMA MCMC workflow for recovered candidates and retain a row for every selected case, including non-detections, wrong aliases and unsupported paths.')
h('Keep independent runs separate')
p('Run four calls with Ncores=1 and distinct recorded seeds rather than immediately concatenating chains. Start with Niter=10000 and retain the cold samples, per-run settings and verbose sampler output through the appropriate mcfit()/sigfit() call. If the data wrapper does not expose a requested sampler control, use the lower-level wrapper explicitly and document the call.')
p('Use trace/rank plots, rank-normalized split Rhat, bulk/tail effective sample sizes and Monte Carlo error [9]. Proposed screening targets are Rhat&lt;1.01 and bulk/tail ESS&gt;400 for each reported scalar. These are diagnostic targets, not a proof of convergence. Inspect multimodality, boundary pile-up and circular parameters separately.')
table(['IMPLEMENTATION DETAIL','CONSEQUENCE FOR THIS TEST'],[
 ('Local period bounds','Single- and multi-set wrappers restrict period to 0.8-1.2 times the starting/reference peak. The prior is flat in sampled log-period, not linear period. Distant aliases cannot be visited.'),
 ('Iterations and burn-in','A positive request is raised to at least 1000 iterations. Default burn-in is floor(iterations/2); extensions can add up to three post-burn-in blocks. Record returned draws, not just Niter.'),
 ('Reported optimum','The current wrappers choose maximum loglike for their point prediction. This is not generally the posterior MAP under a nonuniform eccentricity prior.'),
 ('GP / multiple sets','No single-set GP MCMC. Multi-set GP MCMC holds covariance fixed; multi-set ARMA MCMC does not refine lag coefficients; multi-signal multi-set fits lack a joint posterior refit.'),
], [115,CW-115])
h('Coverage and prior sensitivity')
p('Compute 68% and 95% interval coverage for P, K and e among correctly matched, diagnostically acceptable fits, and report how many cases were excluded. Also report the probability that the whole pipeline both recovers the signal and returns an interval containing truth. Compare Beta(0.867,3.03) and uniform-e priors on a frozen subset.')
p('Fixed-parameter coverage need not equal nominal Bayesian credibility, especially at e=0 or under informative priors. For formal simulation-based calibration, use a separate prior-predictive experiment with fixed, data-independent bounds, matching generator and likelihood, and rank diagnostics [10]. Data-dependent local bounds in the normal wrappers prevent claiming automatic SBC validity.','small')

# 12
page('Stress tests and diagnostics','Additional experiments | Keep these outside the core 3000-case rate denominator')
table(['TEST','CONTROLLED COMPARISON','WHAT TO REPORT'],[
 ('Rotation overlap','Inject periods at 0.5, 0.9, 1.0, 1.1 and 2 times Prot in correlated noise. Pair with un-injected data.','Recovery and amplitude attenuation under W/MA/AR/free GP; compare externally fixed Prot only as a separate informed branch.'),
 ('Eccentric harmonics','Compare one eccentric planet with two circular signals near a 2:1 ratio. Repeat Nh=1 and Nh=2.','Wrong multiplicity, P/2 or 2P confusion, residual peaks and sensitivity to orbital phase coverage.'),
 ('Long-period signal','Use P/T=0.7, 1.2 and 2 with an injected trend; vary seasonal gaps.','Period/trend degeneracy, prediction error and unjustifiably narrow intervals. P/T=2 is an extra boundary case.'),
 ('Instrument changes','Inject known offsets, different jitters and an epoch-dependent step. Fit correct versus deliberately omitted segmentation.','Spurious long periods, offset recovery and the change in combined detection.'),
 ('Sampling loss','Remove whole nights or seasons using a predeclared mask; compare with the complete series.','Loss of sensitivity versus false certainty. Preserve within-night correlation in any bootstrap.'),
 ('Moving diagnostics','Compare a coherent planet, evolving activity and a planet plus activity. Use fixed versus adaptive windows.','Window occupancy, actual widths and diagnostic flags. Do not interpret a normalized moving power as an absolute Bayes factor.'),
], [97,219,CW-316])
p('A manageable supplemental design is 20 paired realizations per stress configuration. Use pilot timings to finalize its size before evaluation. Report counts explicitly; do not silently pool stress cases into the main campaign.')
h('Optional N-body extension')
p('Run tests/test_nbody.R first. Reproduce the code\'s separated stable-pair and close unstable-pair examples, then repeat with half the timestep and a longer short integration. Compare analytic AMD/Hill checks, integration classification and energy error. Start with the R engine; record REBOUND/Python versions if adding that engine.')
note('Keep the claims separate','A stable finite integration does not confirm that an RV signal is planetary or guarantee indefinite stability. Verify that the N-body input comes from the intended accepted 1D fit. Report stellar mass, inclination, timestep, duration and operational instability thresholds.')

# 13
page('Required records and figures','Submission schema | Names are proposed conventions, not existing Agatha2 output fields')
p('Use one row per case and analysis configuration in runs.csv. Use a separate candidates.csv with one row per returned candidate. Keep truth in a separate file and join it only in the scoring stage.')
table(['TABLE','MINIMUM FIELDS'],[
 ('manifest.csv','case_id, partition, input_path, input_hash, template, source, units, N_raw, N_used, N_nights, T, cleaning_version.'),
 ('truth.csv','case_id, seed, generator_version, noise_class, hyperparameters, nplanet, P/K/e/omega/M0 per planet, offsets, trend and jitter. Hidden for validation until scoring.'),
 ('runs.csv','case_id, config_id, commit, settings_hash, seed_fit, status, error_message, warnings_path, seconds, N_requested, N_returned, result_path.'),
 ('candidates.csv','case_id, config_id, candidate_id, P, lnBF, score_source, accepted, matched_truth_id, alias_class, P/K/e estimates and intervals where applicable.'),
 ('diagnostics.csv','case_id, config_id, chain_id, Rhat/ESS/Monte Carlo error, boundary flags, activity/window flags, moving-window counts and failure reason.'),
], [110,CW-110])
h('Figures to include in the final report')
p('<b>1.</b> Sampling/error distributions and the periodograms of the observing windows. <b>2.</b> Completeness versus P/T and K/s, separated by noise and eccentricity. <b>3.</b> Null false-positive rates with intervals and failure fractions. <b>4.</b> True versus fitted P, K and e; alias and multiplicity confusion. <b>5.</b> Interval coverage/width and MCMC diagnostics. <b>6.</b> Paired method differences and runtime scaling. <b>7.</b> One full success and one full failure case per major capability.')
h('Reproducible directory structure')
code('''validation_v1/
  protocol/     # frozen configuration and scoring rules
  inputs/       # originals, cleaned data, manifest, hashes
  truth/        # hidden until scoring
  runs/         # case/config result.rds and logs
  tables/       # run, candidate and score tables
  figures/      # figures plus the data used to draw them
  environment/  # commit, sessionInfo, dependency versions''')
p('Keep generated outputs out of the source commit unless intentionally selected as fixtures. Archive enough inputs and seeds to recreate every result. A CSV alone may lose R attributes such as data-set mappings; preserve the full RDS object too.','small')

# 14
page('Completion and review criteria','Proposed project gates | These assess evidence quality as well as software performance')
table(['GATE','EVIDENCE REQUIRED'],[
 ('Environment ready','The training exercise completes, the strong signal is recovered, and relevant regression failures are resolved or explicitly documented before the campaign.'),
 ('Inputs approved','Every case has provenance, a checksum, validated units/errors and a reproducible cleaning decision. Benchmark scope matches the actual baseline.'),
 ('Protocol frozen','The pilot is separate. Configurations, threshold, selected cases, matching rule and timeout policy are recorded before validation truth is revealed.'),
 ('Jobs accounted for','Every scheduled case/configuration has a result or failure row. Counts reconcile with the matrix; retries preserve the first failure and the changed setting.'),
 ('Statistics correct','Rates state their denominator and uncertainty. Paired methods share the same realization. Aliases and failed jobs are not silently discarded.'),
 ('Claims supported','Conclusions distinguish numerical correctness, sensitivity, false detections, local posterior calibration and observed-target agreement. Unsupported paths are labeled.'),
 ('Independent reproduction','A reviewer reruns five randomly selected cases plus all highlighted failure examples using saved inputs, settings and versions.'),
], [113,CW-113])
h('Suggested final report structure')
p('Use: scientific objective; code/data versions; frozen experiment design; observed benchmarks; null/injection results; posterior calibration; stress tests; implementation limitations; reproducibility instructions. Include a concise capability table with supported, conditional, failed and not tested outcomes.')
h('Scientific performance targets')
p('Assess the proposed 1% null false-positive goal with uncertainty, separately for each major noise regime. Inspect whether detection completeness rises with amplitude in well-covered cases and whether the expected strong benchmark signals are recovered. Investigate violations; do not redefine the sample or threshold after seeing them.')
p('Do not require every weak injection or every planet from the literature to be detected. Do not call biased or undercovered posteriors calibrated because the periodogram peak is correct. If the evidence is insufficient, report the limitation and the additional experiment needed.')
note('A finding is a valid outcome','An honest report of a high false-positive rate, a missed alias, a failed GP fit or a conditional multi-set uncertainty is a useful validation result. The assessment is complete when the result is reproducible and its scope is stated, not when every plot looks favorable.')

# 15
page('Troubleshooting and starting results','Executed training exercise | Core code snapshot 229dccc')
table(['SYMPTOM','STUDENT RESPONSE'],[
 ('Missing functions/packages','Confirm repository working directory, package installation and source order. Preserve the error and sessionInfo before changing the environment.'),
 ('No accepted signal','Inspect model.comp before the plotted fit. Determine whether the case is below threshold, an alias, an optimizer failure or an unsupported score.'),
 ('Very slow GP or MCMC','Use pilot timings, queued work and an explicit core budget. Reduce a balanced design before the final freeze; do not selectively drop difficult validation cases.'),
 ('Posterior misses another peak','Check local period bounds and the initialization. A distinct starting peak is a separate local analysis; mixing separate chains without mode weights is not a global posterior.'),
 ('MCMC absent','Check Niter and model support. Single-set GP MCMC is unavailable; stochastic fits do not return an orbital posterior.'),
 ('Results change across runs','Check code/data hashes, seeds, parallel RNG behavior, package versions and any implicit loader transformations. Repeat with Ncores=1 for diagnosis.'),
 ('Unexpected extra planet','Inspect harmonics, activity, window aliases, trends and instrument steps. Retain the result in the score tables.'),
], [119,CW-119])
h('Training execution during report preparation')
p('The companion script was run on this checkout on 03 October 2026. It used 80 synthetic observations over approximately 180 days, unit quoted uncertainties and a 5 m/s circular injection at 17.3 days. The null and injected series shared the same white-noise realization.')
table(['CASE','RETURNED PEAK','ln(BF)','ACCEPTED'],[
 ('Noise only','10.418 days','0.226','No'),
 ('Injected signal','17.268 days','98.315','Yes'),
], [120,155,104,CW-379])
p('These are actual training outputs, not predictions for the full campaign. Their elapsed times were about 32 and 18 seconds on the preparation machine; use new local pilot timings for resource planning. The output directory includes the synthetic inputs and the exact settings.','small')
note('Work still to be done','Students must implement and verify the full generator, scheduler and scorer; run the independent null/injection campaign; audit the challenge truth; and obtain any optional external benchmark data. No full scientific validation or new long MCMC/N-body campaign was completed while writing this report.')

# 16
page('References and reproducibility','Primary sources and local implementation | Prepared 03 October 2026')
p('Use the original papers and data releases to define benchmark expectations. The protocol\'s allocation, thresholds and diagnostic gates are proposed experimental choices, not claims made by those papers.','small')
refs=[
 ('1','Butler et al. (2004). A Neptune-Mass Planet Orbiting the Nearby M Dwarf GJ 436.','https://arxiv.org/abs/astro-ph/0408587'),
 ('2','Mayor et al. (2009). The HARPS search for southern extra-solar planets. XIII. A planetary system with 3 Super-Earths.','https://arxiv.org/abs/0806.4587'),
 ('3','Unger et al. (2021). The HARPS search for southern extra-solar planets. XLVI. Includes HD 189567.','https://research-repository.st-andrews.ac.uk/handle/10023/24291'),
 ('4','Feng, Jones &amp; Tuomi (2018). An improved quantification of HD 147379 b.','https://arxiv.org/abs/1801.10541'),
 ('5','Tamuz et al. (2008). Discovery of two eccentric planets orbiting HD 4113 and HD 156846.','https://arxiv.org/abs/0710.5028'),
 ('6','Queloz et al. (2001). No planet for HD 166435.','https://arxiv.org/abs/astro-ph/0109491'),
 ('7','Dumusque (2016). Radial Velocity Fitting Challenge. I. Data-generation and injected-system definitions.','https://arxiv.org/abs/1607.06487'),
 ('8','Dumusque et al. (2021). Three years of HARPS-N high-resolution spectroscopy and precise radial velocity data for the Sun.','https://www.aanda.org/articles/aa/full_html/2021/04/aa39350-20/aa39350-20.html'),
 ('9','Vehtari et al. (2021). Rank-normalization, folding, and localization: An improved Rhat for assessing convergence of MCMC.','https://arxiv.org/abs/1903.08008'),
 ('10','Talts et al. (2018; revised 2020). Validating Bayesian Inference Algorithms with Simulation-Based Calibration.','https://arxiv.org/abs/1804.06788'),
]
for n,title,url in refs:
    p(f'[{n}] {title} <link href="{url}" color="#007F86">Source</link>.','small')
h('Local source of implementation claims')
p('functions.R: calc.1Dper(), mcfit(), mcfit.multiset(), read.agatha.table(), diagnose.signals(). periodoframe.R: BFP(), GP and multi-set fitting. mcmc_func.R: run.ptmcmc(). additional_signals.R: sequential acceptance. nbody.R: stability utilities. data/DACE_TEST_DATA.md: archived input provenance. tests/: existing regression checks.','small')
h('Companion files and document build')
code('''Rscript docs/examples/validation_smoke.R
python3 docs/build_student_validation.py''')
p('The PDF builder requires Python and reportlab. The starter requires the R dependencies on page 3. Distribute the exact source checkout, the starter and this protocol together. Agatha2 is distributed under its MIT license; preserve license.txt with redistributed source.','small')

doc=GuideDoc(str(OUT),pagesize=A4,leftMargin=M,rightMargin=M,
    topMargin=58,bottomMargin=58,title='Agatha2 - Student Validation Protocol',
    author='Agatha2 project',subject='Student experiments for detection and uncertainty validation')
doc.addPageTemplates([PageTemplate(id='guide',frames=[Frame(M,58,CW,H-116,id='body',
    leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=footer)])
doc.build(story)
print(OUT)
