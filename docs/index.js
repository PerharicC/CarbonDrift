// Load docs

const navBarWidthThreshold = 768;

document.addEventListener("DOMContentLoaded", function () {
    if (document.title != "CarbonDrift") {
        if (window.innerWidth < navBarWidthThreshold) {
            $(".navbar").attr("hidden", true);
        } else {
            $(".navbar-collapse").attr("hidden", true);
            $("#btn-collapse-navbar").attr("hidden", true)
        };
        fetch("./navbar.html")
            .then(response => response.text())
            .then(data => {
                document.querySelector(".navbar").innerHTML = data;
            })
        fetch("./collapsednavbar.html")
            .then(response => response.text())
            .then(data => {
                document.querySelector(".navbar-collapse").innerHTML = data;
            })
        setTimeout(()=>{buttonChange()}, 100);
        
        const plotId = localStorage.getItem("plotId");
        if (plotId) {
            expandMethod($(`#${plotId}`)[0]);
            localStorage.removeItem(plotId);
        };
    } else {
        let cards = $(".card");
        let maxWidth = 0
        for (let i = 0; i < cards.length; i++) {
            maxWidth = Math.max($(cards[i]).width(), maxWidth);
        };
        setTimeout(()=>{
            $(".card").css("width", maxWidth);
            console.log($(".card").width())
        }, 100);
    };
});

// Get Current Page in Docs
function getPage() {
    switch (document.title) {
        case "Theory":
            newId = "btn-theory";
            break;
        case "Install":
            newId = "btn-installation";
            break;
        case "Terminal Simulation":
            newId = "btn-terminal";
            break;
        case "Custom Simulation":
            newId = "btn-custom-simulation";
            break;
        case "GUI":
            newId = "btn-gui";
            break;
        case "Blueprint":
            newId = "btn-blueprint";
            break;
        case "Plot":
            newId = "btn-plotting";
            break;
        case "Gallery":
            newId = "btn-gallery";
            break;
        case "Release Notes":
            newId = "btn-release-notes";
            break;
        case "API References":
            newId = "btn-api";
            break;
        case "References":
            newId = "btn-references";
            break;
        default:
            newId = "btn-theory"
    };
    return newId;
};

// Update NavBar
function buttonChange(){
    newId = getPage();
    document.querySelectorAll(".nav-link").forEach(link => {
                    link.classList.remove("active");
                    link.classList.add("link-body-emphasis")
                    link.removeAttribute("aria-current");
                });
    
    let newHighlitedButton = document.getElementById(newId);
    newHighlitedButton.classList.add("active");
    newHighlitedButton.classList.remove("link-body-emphasis");
    newHighlitedButton.setAttribute("aria-current", "page");

    let newHighlitedButton2 = document.getElementById(newId+"2");
    newHighlitedButton2.classList.add("active");
    newHighlitedButton2.classList.remove("link-body-emphasis");
    newHighlitedButton2.setAttribute("aria-current", "page");
}

$(window).on("resize", (e) =>{
    let windowWidth = e.target.innerWidth;
    if (windowWidth>=navBarWidthThreshold) {
        $(".navbar-collapse").attr("hidden", true);
        $("#btn-collapse-navbar").attr("hidden", true)
        $(".navbar").removeAttr("hidden");
    } else {
        $(".navbar").attr("hidden", true);
        $(".navbar-collapse").removeAttr("hidden");
        $("#btn-collapse-navbar").removeAttr("hidden");
    };
});

//Gallery

$(".gallery-button").mouseover(function (e) {
    e.target.classList.toggle("hover")
});

$(".gallery-button").mouseout(function (e) {
    e.target.classList.remove("hover");
});

$("#gallery img").click(function(e) {
    largeImage = e.target;
    let gifSrc = e.target.src + "?t=" + new Date().getTime();
    $("#lightbox-img")[0].src = gifSrc;
    $("#lightbox")[0].style.display = "flex";
    $("#lightbox").removeAttr("hidden");
});

$(".gallery-button").click(function (e) {
    let plotId = e.target.id.replace("btn-", "").replace("-copy", "")
    localStorage.setItem("plotId", plotId)
    $(".expanded-code .code-method").attr("href", "apiplots.html#" + plotId);
    $(".expanded-code .code-method")[0].innerHTML = getMethodName(plotId);
    let card = $(e.target).parents().eq(2);
    $("#expanded-img")[0].src = card.find("img")[0].src;
    $("#gallery > .container").attr("hidden", true)
    $(".expanded-code h2")[0].innerHTML = card.find(".card-text")[0].innerHTML;
    $(".expanded-code pre")[0].innerHTML = getParameterFile(e.target.id.replace("btn-", ""))
    $("#expanded-code-caption")[0].innerHTML = `This is an example showing how to create a ${getMethodDescription(e.target.id.replace("btn-", ""))}`
    $(".expanded-code")[0].style.display = "block"
})

$(".expanded-code button").click(function () {
    $("#expanded-img")[0].src = ""
    $("#gallery > .container").removeAttr("hidden")
    $(".expanded-code")[0].style.display = "none"
})

function getMethodName(id) {
    switch (id) {
        case "mass-flux-map":
            methodName = "mass_flux_map";
            break;
        case "drifter-locations":
            methodName = "drifter_locations";
            break;
        case "drifter-properties":
            methodName = "drifter_properties";
            break;
        case "flux-distribution":
            methodName = "mass_flux_distribution";
            break;
        case "lat-mean":
            methodName = "mean_lon_mass_flux";
            break;
        case "biome-flux-mean":
            methodName = "mean_export_biome_flux";
            break;
        case "current-strength":
            methodName = "current_strength";
            break;
        case "sinking-animation":
            methodName = "animate_3D";
            break;
        case "current-animation":
            methodName = "animate_current_3D";
            break;
        case "mass-table":
            methodName = "get_biome_weighted_mass_at_depth";
            break;
        case "z-dist":
            methodName = "vertical_particle_distribution";
            break;
        default:
            methodName = "";
    };
    return methodName;
}

function getParameterFile(parameterID){
    switch (parameterID) {
        case "mass-flux-map":
            text = "\tmass_flux_map\n\tfile1\n\tfile2\n\t...\n\tfileN\n\t-fs\n\t12:8\n\t-cblabel\
            \n\t$\\text{carbon flux at 100 m } \\mathrm{[g\\,C\\,m^{-2}\\,Y^{-1}]}$\n\t-t\n\tEXP_decay,constant_speed\
            \n\t-d\n\t-6000\n\t-cmap\n\ttab20b\n\t--shrink\n\t0.8\n\t-add\
            \n\t-clip\n\t0:20\n\t-fts\n\t25";
            break;
        case "mass-flux-map-copy":
            text = "\tmass_flux_map\n\tfile1\n\tfile2\n\t...\n\tfileN\n\t-fs\n\t12:8\n\t-cblabel\
            \n\t$\\text{Seafloor difference } \\mathrm{[g\\,C\\,m^{-2}\\,Y^{-1}]}$\n\t-t\n\tSSP585\
            \n\t-d\n\t-6000\n\t-cmap\n\ttab20b\n\t--shrink\n\t0.8\n\t-add\n\t-diff\n\t-abs\
            \n\t-clip\n\t1:3\n\t-fts\n\t25";
            break;
        case "drifter-locations":
            text = "\tdrifter_locations\n\tfile1\n\t-fs\n\t20:10\n\t-loc\
            \n\tm48:30,m2:75,1:38\n\t-legend\n\tMid Atlantic,Arctic,Western Mediterranean\
            \n\t-color\n\tblue,black,red\n\t-xlim\n\tm100:50\n\t-ylim\n\t0:90\n\t-fts\n\t22";
            break;
        case "drifter-properties":
            text = "\tdrifter_properties\n\tfile1\n\tfile2\n\t...\n\tfileN\n\t-fs\n\t20:10\n\t-loc\
            \n\tm48:30,m2:75,1:38\n\t-p1\n\tmass\n\t-p2\n\tz\n\t-xlabel\n\t$m/m_0$\n\t-ylabel\n\t$z\,\\mathrm{[m]}$\n\t-legend\
            \n\tEXP $\\frac{dm}{dt} = -km$,LIN $\\frac{dm}{dt} = -km$,EXP $\\frac{dm}{dt} = -kS$,LIN $\\frac{dm}{dt} = -kS$\
            \n\t-color\n\troyalblue,blue,darkorange,orange\n\t-xlim\n\t0:1\n\t-ylim\n\tm2500:0\n\t-lw\n\t3\n\t-ls\n\tsolid,dashed,solid,dashed\
            \n\t-mc\n\t0.86\n\t-supt\n\tMid Atlantic,Arctic,Western Mediterranean\n\t-fts\n\t19";
            break;
        case "flux-distribution":
            text = "\tmass_flux_distribution\n\tfile1\n\tfile2\n\t...\n\tfileN\n\t-fs\n\t8:10\
            \n\t-add\n\t-group\n\tbiome\n\t-d\n\t-100\n\t-xlabel\n\t$\\text{mass flux at 100 m}\\mathrm{[g\\,C\\,Y^{-1}\\,m^{-2}]}$\n\t-t\
            \n\tEXP_Area_Decay\
            \n\t-color\n\tblue,red,orange,green\n\t-lw\n\t2\n\t-bins\n\t10\
            \n\t-xlim\n\t0:30\n\t-fts\n\t20";
            break;
        case "lat-mean":
            text = "\tmean_lon_mass_flux\n\tfile1\n\tfile2\n\t...\n\tfileN\n\t-fs\n\t6:6\
            \n\t-add\n\t-legend\n\tmass EXP,mass LIN,area EXP,area LIN\n\t-d\n\t-6000\n\t-ylabel\n\tLatitude averaged seafloor flux $\\mathrm{[g\\,C\\,m^{-2}\\,Y^{-1}]}$\n\t-xlabel\
            \n\t$\\varphi$\
            \n\t-color\n\tblue,red,blue,red\n\t-lw\n\t2.5\n\t-ls\n\tsolid,solid,dashed,dashed\
            \n\t-ylim\n\t0:25\n\t-fts\n\t16";
            break;
        case "biome-flux-mean":
                text = "\tmean_export_biome_flux\n\tfile1\n\tfile2\n\t...\n\tfileN\n\t-fs\n\t20:6\
                \n\t-legend\n\tChordata,Ctenophora,Cnidaria\n\t-d\n\t-100\n\t-ylabel\n\tmass flux at 100 m $\\mathrm{[g\\,C\\,m^{-2}\\,Y^{-1}]}$\n\t-xlabel\
                \n\tEXP constant,EXP variable,LIN constant,LIN variable\
                \n\t-color\n\tblue,orange,magenta\n\t-lw\n\t2.5\n\t-ls\n\tsolid,solid,dashed,dashed\
                \n\t-supt\n\tmortality,egestion";
                break;
        case "current-strength":
            text = "\tcurrent_strength\n\tfile1\n\t-fs\n\t12:8\
            \n\t-cblabel\n\t$|\\vec{r}_H(z=-1000\\,\\mathrm{m}) - \\vec{r}_{H0}|$ [km]\n\t-d\n\t-1000\n\t--shrink\
            \n\t0.8\n\t-cmap\n\tBlues\n\t-t\n\tw_0_=_100_m/d\
            \n\t-clip\n\t0:100\n\t-fts\n\t30";
            break;
        case "sinking-animation":
            text = "\tanimate_3D\n\tfile1\n\t-fs\n\t12:8\
            \n\t--shrink\
            \n\t0.6\n\t-fts\n\t21";
            break;
        case "sinking-animation-copy":
            text = "\tanimate_3D\n\tfile1\n\t-fs\n\t12:8\
            \n\t-cmap\n\tmagma\n\t--shrink\
            \n\t0.6\n\t-fts\n\t21";
            break;
        case "current-animation":
            text = "\tanimate_current_3D\n\tfile1\n\t-fs\n\t15:8\
            \n\t-xlim\n\tm75:m62\n\t-ylim\n\tm67:m53\n\t-cmap\n\tcool\n\t--shrink\
            \n\t0.6\n\t-fts\n\t21";
            break;
        case "mass-table":
            text = "\tget_biome_weighted_mass_at_depth\n\tfile1\n\t-abs\n\t-d\n\t-100";
            break;
        case "z-dist":
            text = "\tvertical_particle_distribution\n\tfile1\n\t...\n\tfileN\n\t-fs\n\t8:10\n\t-color\
            \n\tred,orange,lightcoral,navy,royalblue,skyblue\n\t-bins\n\t45\n\t-t\n\tParticle Distribution\
            \n\t-d\n\t-4500\n\t-legend\n\tCnidaria $E_g$,Ctenophora $E_g$,Chordata $E_g$,Cnidaria $M$,Ctenophora $M$,Chordata $M$\
            \n\t-xlim\n\t0:0.25\n\t-xlabel\n\trelative number of particles\n\t-ylabel\n\tz [m]";
            break;
        default:
            text = "";
    };
    return text;
};

function getMethodDescription(parameterID){
    switch (parameterID) {
        case "mass-flux-map":
            text = "mass flux map at a given depth by stacking multiple simulations together.";
            break;
        case "mass-flux-map-copy":
            text = "difference map between mass fluxex at a given depth. The first N//2 and the rest of the simulations\
            are stacked together, respectively. Next, their differnce is computed.";
            break;
        case "drifter-locations":
            text = "plot of speciffic locations.";
            break;
        case "drifter-properties":
            text = "plot of different dynamical properties at specified locations and for different simulations. \
            One can choose between all variables, which are saved to the netCDF file. Additionaly, one can compute the the \
            total horizontal velocity norm or the inverse decay rate by setting -p1/-p2 to total_horizontal_velocity or tau, respectively.";
            break;
        case "flux-distribution":
            text = "distribution plot of mass fluxes in each biome.";
            break;
        case "lat-mean":
            text = "mass flux averaged over all longitudes at a given latitude plot.";
            break;
        case "biome-flux-mean":
            text = "average flux in each biome for multiple simulations. The order of the files should be \
            exp_variable_chordata_M, exp_variable_chordata_Eg, exp_variable_ctenophora_M, exp_variable_ctenophora_Eg, exp_variable_cnidaria_M, exp_variable_cnidaria_Eg\
            and so on for the next type of decay...";
            break;
        case "current-strength":
            text = "plot of horizontal distances traveled by tracers from their initial positions to their positions at a given depth.";
            break;
        case "sinking-animation":
            text = "sinking animation.";
            break;
        case "sinking-animation-copy":
            text = "sinking animation.";
            break;
        case "current-animation":
            text = "sinking and current strength animation.";
            break;
        case "mass-table":
            text = "tabulated value of the total mass of particles at a given depth. In order to retrieve a full table as shown in the \
            image example, run the following script for multiple depths and for different simulations (phylums).";
            break;
        case "z-dist":
            text = "vertical particle distribution animation. If the outfile is not given, an interactive slider animation will be created.";
            break;
        default:
            text = "";
    };
    return text;
};

//API

function expandMethod(target) {
    if (typeof target !== "undefined") {
        const content = target.querySelector("div")
        content.classList.toggle("method-content")
        target.classList.toggle("expanded")
    };
}

$(".expandable > li").click(function (e){
    expandMethod(e.target);
});

$(".gallery-button").click(function (e) {
    localStorage.setItem("plotId", e.target.id.replace("btn-", "").replace("-copy", ""));
});

// Footer

function copyrightText() {
    const currentDate = new Date();
    let currentYear = currentDate.getFullYear();
    return `Copyright ${currentYear}, Črtomir E. Perharič Bailey`
}

document.addEventListener("DOMContentLoaded", function () {
    setTimeout(()=>{
        $("#copyright-text")[0].innerHTML = copyrightText()
    }, 100);
});

document.addEventListener("DOMContentLoaded", function () {
    if (document.title != "CarbonDrift") {
        fetch("./footer.html")
            .then(response => response.text())
            .then(data => {
                document.querySelector(".footer").innerHTML = data;
            });
        setTimeout(()=>{
            $("#copyright-text")[0].innerHTML = copyrightText()
        }, 100);
    };
    });

// Main Page

$("#CD-summary p").css("width", `${window.innerWidth * 0.007 +493}px`)

$(window).on("resize", () => {
    if (document.title === "CarbonDrift") {
        let cards = $(".card");
        let maxWidth = 0
        for (let i = 0; i < cards.length; i++) {
            maxWidth = Math.max($(cards[i]).width(), maxWidth);
        };
        $(".card").css("width", maxWidth);
    }
})

$(".btn-header").mouseover(function (e) {
    e.target.classList.toggle("btn-header-hover");
});

$(".btn-header").mouseout(function (e) {
    e.target.classList.remove("btn-header-hover");
});

$("#btn-header-features").click(function(e) {
    $("#btn-header-home")[0].classList.remove("btn-header-active");
    $("#btn-header-features")[0].classList.add("btn-header-active");
    $("#image-examples").attr("hidden", true);
    $("#features").removeAttr("hidden");
    let search = document.getElementById("searched-elements");
    if (! search.hasAttribute("hidden")) {
        toggleSections(false)
        search.setAttribute("hidden", true)
    };
});

$("#btn-header-home").click(function(e) {
    $("#btn-header-features")[0].classList.remove("btn-header-active");
    $("#btn-header-home")[0].classList.add("btn-header-active");
    $("#features").attr("hidden", true);
    $("#image-examples").removeAttr("hidden");
    let search = document.getElementById("searched-elements");
    if (! search.hasAttribute("hidden")) {
        toggleSections(false)
        search.setAttribute("hidden", true)
    };
});

$("#btn-header-docs").click(function (){
    localStorage.setItem("newLink", "btn-theory");
})

$("#image-examples img").click(function(e) {
    largeImage = e.target;
    let gifSrc = e.target.src + "?t=" + new Date().getTime();
    $("#lightbox-img")[0].src = gifSrc;
    $("#lightbox")[0].style.display = "flex";
    $("#lightbox").removeAttr("hidden");
});

$(".lightbox-close").click(function() {
    $("#lightbox")[0].style.display = "none";
    $("#lightbox").attr("hidden", true);
});

// Search

$("#search-bar").on("submit", async function(e){
    e.preventDefault();
    const searchInput = $("input").first().val().toLowerCase();
    let matchedSearchElements = await searchElements(searchInput);
    $("#search input")[0].value = "";
    if (matchedSearchElements.length === 0) {
        $("#search input").attr("placeholder", "No match. Try again...");
    } else {
        displaySearchedElements(matchedSearchElements);
    };
});

async function searchElements(searchInput) {
    try {
        const response = await fetch("./search.json");
        const data = await response.json();
        let matches = data.filter(item => (
                item.title.split(" ").filter(opt => (
                    opt.includes(searchInput) || searchInput.includes(opt)
                ))
            ).length >0);
        return matches;
        } catch (error) {
            console.error('Error loading search data:', error);
            return [];
    };
}

function toggleSections(hide) {
    let sections = $("section");
    for (let i = 0; i < sections.length; i++) {
        if (! sections[i].classList.contains("title")) {
            let hidden = $(sections[i]).attr("hidden");
            if (typeof hidden === "undefined" || hidden === false || hide === true) {
                $(sections[i]).attr("hidden", true);
            } else {
                $(sections[i]).removeAttr("hidden");
            };
        };
      };
}

function displaySearchedElements(elements) {
    $("#searched-elements ul").remove();
    $("#searched-elements").append($("<ul></ul>"))
    elements.forEach(item => {
        let searchLink = $("<a></a>").attr("href", item.url).text(item.display);
        let searchSVG = $("<img />").attr("src", "./assets/svgs/search.svg").addClass("search-list")
        let searchListElement = $("<li></li>")
        searchListElement.append(searchSVG);
        searchListElement.append(searchLink);
        $("#searched-elements ul").append(searchListElement);
    })
    toggleSections(true);
    $("#searched-elements").removeAttr("hidden")
}

// Copy Code Button

let codeBlocks = $("pre")
let copySVG = "./assets/svgs/copy.svg"
codeBlocks.wrap("<div class='code-div'></div>")
if (navigator.clipboard) {
    let button = $("<button></button>")
    button.append($(`<img src= "${copySVG}" />`))
    button.wrap("<div class='copy-code'></div>")
    $("pre").before(button.parent())
}

$(".copy-code button").click((e)=> {
    let codeWrap = e.target.closest(".code-div")
    let codeBlock = $(codeWrap).find("code")
    if (codeBlock.length === 0) {
        codeBlock = $(codeWrap).find("pre")
    }
    let code = codeBlock.text()
    navigator.clipboard.writeText(code)
})
