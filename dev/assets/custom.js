// document.addEventListener("DOMContentLoaded", function () {
//     localStorage.setItem("documenter-theme", "catppuccin-mocha");
// });

(function () {
    const theme = "documenter-dark";

    // Apply the theme before the page renders
    if (localStorage.getItem("documenter-theme") !== theme) {
        localStorage.setItem("documenter-theme", theme);
    }

    // Enforce the theme immediately on page load
    document.documentElement.setAttribute("data-theme", theme);
})();

// document.addEventListener("DOMContentLoaded", function () {
//     document.querySelectorAll("html.theme--catppuccin-mocha a").forEach(el => {
//         el.addEventListener("mouseover", function () {
//             if (getComputedStyle(this).color === "rgb(137, 220, 235)") { // #89dceb in RGB
//                 this.style.color = "#4493f8";
//             }
//         });

//         el.addEventListener("mouseout", function () {
//             this.style.color = ""; // Resets to default when not hovered
//         });
//     });
// });