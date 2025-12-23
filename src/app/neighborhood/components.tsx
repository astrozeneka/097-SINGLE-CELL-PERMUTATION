

function Form() {
  return (
    <div>
      <h2>Neighborhood Analysis</h2>
      <div>
        <label htmlFor="method-select">Select Method:</label>
        <select id="method-select">
            <option value="knn">k-Nearest Neighbors</option>
            <option value="radius">Radius-based</option>
        </select>
      </div>
      <div>
        {/* In case of knn, show k input */}
        <div>
          <label htmlFor="k-input">Number of Neighbors (k):</label>
          <input type="number" id="k-input" min="1" defaultValue="5" />
        </div>
        {/* In case of radius, show radius input */}
        <div>
          <label htmlFor="radius-input">Radius:</label>
          <input type="number" id="radius-input" min="0.1" step="0.1" defaultValue="10.0" />
        </div>
      </div>
    </div>
  );
}